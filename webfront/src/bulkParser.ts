export const NOSTRUCT_COLUMNS = [
  "External ID",
  "External Batch",
  "MW",
  "restriction_comment"
] as const;

export interface SdfImport {
  kind: "sdf";
  columns: string[];
  rows: Array<{ molfile: string; fields: Record<string, string> }>;
}

export interface SmilesImport {
  kind: "smiles";
  columns: string[];
  rows: Array<Record<string, string>>;
}

export interface NostructImport {
  kind: "nostruct";
  rows: Array<{
    externalId: string;
    supplierBatch: string;
    mw: number;
    restrictionComment: string;
  }>;
}

export type ParsedBulkFile = SdfImport | SmilesImport | NostructImport;

function cellText(value: unknown): string {
  if (value === null || value === undefined) return "";
  if (typeof value === "string" || typeof value === "number" || typeof value === "boolean") {
    return String(value).trim();
  }
  if (typeof value === "object") {
    const cell = value as { text?: unknown; richText?: Array<{ text?: unknown }> };
    if (typeof cell.text === "string") return cell.text.trim();
    if (Array.isArray(cell.richText)) {
      return cell.richText.map((item) => String(item.text ?? "")).join("").trim();
    }
  }
  return "";
}

function parseSdfRecord(block: string): SdfImport["rows"][number] {
  const lines = block.replace(/\r\n/g, "\n").replace(/\n+$/, "").split("\n");
  if (!lines[0]?.trim()) lines[0] = "id";
  const fields: Record<string, string> = {};

  for (let index = 0; index < lines.length; index += 1) {
    const match = /^>\s*<([^>]+)>/.exec(lines[index] ?? "");
    if (!match) continue;
    const values: string[] = [];
    index += 1;
    while (index < lines.length && lines[index]?.trim()) {
      values.push(lines[index] ?? "");
      index += 1;
    }
    fields[match[1].trim()] = values.join("\n").trim();
  }

  return { molfile: `${lines.join("\n")}\n$$$$\n`, fields };
}

export function parseSdfText(text: string): SdfImport {
  const blocks = text
    .replace(/\r\n/g, "\n")
    .split(/^\$\$\$\$\s*$/m)
    .filter((block) => block.trim());
  const rows = blocks.map(parseSdfRecord);
  if (!rows.length) throw new Error("The SDF file contains no molecule records.");
  const columns = [...new Set(rows.flatMap((row) => Object.keys(row.fields)))];
  return { kind: "sdf", columns, rows };
}

export function validateNostructRows(
  rows: NostructImport["rows"]
): string[] {
  const errors: string[] = [];
  if (!rows.length) return ["The workbook must contain at least one data row."];

  rows.forEach((row, index) => {
    const rowNumber = index + 2;
    if (!row.externalId) errors.push(`Row ${rowNumber}: External ID is required.`);
    else if (row.externalId.length > 200) errors.push(`Row ${rowNumber}: External ID must be at most 200 characters.`);
    if (!row.supplierBatch) errors.push(`Row ${rowNumber}: External Batch is required.`);
    else if (row.supplierBatch.length > 100) errors.push(`Row ${rowNumber}: External Batch must be at most 100 characters.`);
    if (!Number.isFinite(row.mw)) errors.push(`Row ${rowNumber}: MW must be a finite number.`);
  });
  return errors;
}

async function parseWorkbook(buffer: ArrayBuffer): Promise<ParsedBulkFile> {
  const XLSX = await import("@e965/xlsx");
  const workbook = XLSX.read(buffer, { type: "array", cellDates: false });
  const firstSheetName = workbook.SheetNames[0];
  if (!firstSheetName) throw new Error("The workbook has no worksheets.");
  const sheet = workbook.Sheets[firstSheetName];
  if (!sheet) throw new Error("The first worksheet could not be read.");

  const grid = XLSX.utils.sheet_to_json<unknown[]>(sheet, {
    header: 1,
    raw: true,
    defval: "",
    blankrows: false
  });
  const columns = (grid[0] ?? []).map(cellText);
  while (columns.length && !columns[columns.length - 1]) columns.pop();
  const dataRows = grid.slice(1).filter((row) => row.some((value) => cellText(value)));

  if (columns.length === NOSTRUCT_COLUMNS.length &&
      NOSTRUCT_COLUMNS.every((column, index) => columns[index] === column)) {
    const rows = dataRows.map((row) => {
      const mwText = cellText(row[2]);
      return {
        externalId: cellText(row[0]),
        supplierBatch: cellText(row[1]),
        mw: mwText ? Number(mwText) : 0,
        restrictionComment: cellText(row[3])
      };
    });
    const errors = validateNostructRows(rows);
    if (errors.length) throw new Error(errors.join("\n"));
    return { kind: "nostruct", rows };
  }

  if (!columns.includes("Smile")) {
    throw new Error('The workbook must contain a "Smile" column or the exact no-structure columns.');
  }
  if (!dataRows.length) throw new Error("The workbook must contain at least one data row.");

  const rows = dataRows.map((row) => Object.fromEntries(
    columns.map((column, index) => [column, cellText(row[index])])
  ));
  return { kind: "smiles", columns, rows };
}

export async function parseBulkFile(file: File): Promise<ParsedBulkFile> {
  const extension = file.name.split(".").pop()?.toLowerCase();
  if (extension === "sdf") return parseSdfText(await file.text());
  if (extension === "xlsx" || extension === "xls") {
    return parseWorkbook(await file.arrayBuffer());
  }
  throw new Error("Choose an SDF or Excel workbook.");
}