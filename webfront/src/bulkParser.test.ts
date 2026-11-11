import * as XLSX from "@e965/xlsx";
import { describe, expect, it } from "vitest";
import { parseBulkFile, parseSdfText } from "./bulkParser";

describe("bulk import parsing", () => {
  it("splits SDF records and extracts the first value for property tags", () => {
    const result = parseSdfText([
      "molecule one",
      "  ChemReg",
      "  0  0  0     0  0            999 V2000",
      "M  END",
      ">  <External ID>",
      "EXT-1",
      "",
      "$$$$",
      "molecule two",
      "  ChemReg",
      "  0  0  0     0  0            999 V2000",
      "M  END",
      ">  <External ID>",
      "EXT-2",
      "",
      "$$$$"
    ].join("\n"));

    expect(result.kind).toBe("sdf");
    expect(result.rows).toHaveLength(2);
    expect(result.columns).toEqual(["External ID"]);
    expect(result.rows[0]?.fields["External ID"]).toBe("EXT-1");
    expect(result.rows[0]?.molfile).toContain("$$$$");
  });

  it("recognizes Smile workbooks and maps values to column names", async () => {
    const workbook = XLSX.utils.book_new();
    XLSX.utils.book_append_sheet(
      workbook,
      XLSX.utils.aoa_to_sheet([["Smile", "Supplier Ref"], ["CCO", "A-12"]]),
      "Compounds"
    );
    const data = XLSX.write(workbook, { type: "array", bookType: "xlsx" });
    const file = new File([data], "smiles.xlsx");
    const parsed = await parseBulkFile(file);

    expect(parsed.kind).toBe("smiles");
    if (parsed.kind !== "smiles") return;
    expect(parsed.columns).toEqual(["Smile", "Supplier Ref"]);
    expect(parsed.rows[0]).toEqual({ Smile: "CCO", "Supplier Ref": "A-12" });
  });

  it("reads legacy XLS workbooks used by the desktop SMILES workflow", async () => {
    const workbook = XLSX.utils.book_new();
    XLSX.utils.book_append_sheet(
      workbook,
      XLSX.utils.aoa_to_sheet([["Smile", "External ID"], ["CCO", "EXT-1"]]),
      "Compounds"
    );
    const data = XLSX.write(workbook, { type: "array", bookType: "xls" });
    const parsed = await parseBulkFile(new File([data], "legacy.xls"));

    expect(parsed.kind).toBe("smiles");
    if (parsed.kind !== "smiles") return;
    expect(parsed.rows[0]).toEqual({ Smile: "CCO", "External ID": "EXT-1" });
  });

  it("validates exact no-structure worksheets and treats blank MW as zero", async () => {
    const workbook = XLSX.utils.book_new();
    XLSX.utils.book_append_sheet(
      workbook,
      XLSX.utils.aoa_to_sheet([
        ["External ID", "External Batch", "MW", "restriction_comment"],
        ["EXT-1", "BATCH-1", "", "No structure"]
      ]),
      "Compounds"
    );
    const file = new File([XLSX.write(workbook, { type: "array", bookType: "xlsx" })], "controls.xlsx");
    const parsed = await parseBulkFile(file);

    expect(parsed.kind).toBe("nostruct");
    if (parsed.kind !== "nostruct") return;
    expect(parsed.rows[0]?.mw).toBe(0);
    expect(parsed.rows[0]?.externalId).toBe("EXT-1");
  });
});