import { useEffect, useState, type ChangeEvent } from "react";
import {
  AlertCircle,
  Check,
  Download,
  FileUp,
  LoaderCircle,
  TriangleAlert,
  Upload
} from "lucide-react";
import {
  ApiError,
  getColumnValues,
  getLastBatchOfEln,
  getLastCtrlCompound,
  getNextSdfSequence,
  getRegnosFromSdfSequence,
  postChemRegMol,
  postCtrlMol,
  postNostructMol,
  registerCompound,
  type Session
} from "./api";
import { parseBulkFile, type ParsedBulkFile } from "./bulkParser";

interface UploadMetadata {
  chemist: string;
  compoundType: string;
  project: string;
  supplier: string;
  solvent: string;
  product: string;
  library: string;
  ipRights: string;
}

interface ImportIssue {
  row: string;
  identifier: string;
  kind: "error" | "warning";
  message: string;
}

interface ImportReport {
  total: number;
  submitted: number;
  registered: number;
  newCompounds: number;
  issues: ImportIssue[];
}

interface BulkImportScreenProps {
  session: Session;
  onSessionExpired: () => void;
}

const emptyMetadata: UploadMetadata = {
  chemist: "",
  compoundType: "",
  project: "",
  supplier: "",
  solvent: "",
  product: "",
  library: "",
  ipRights: ""
};

const ipRightsOptions = ["External rights", "LCBKI", "Commercial"];

function errorText(error: unknown): string {
  return error instanceof Error ? error.message : "The request could not be completed.";
}

function chooseDefaultColumn(columns: string[], matches: string[]): string {
  const lowered = new Map(columns.map((column) => [column.toLowerCase(), column]));
  for (const match of matches) {
    const column = lowered.get(match.toLowerCase());
    if (column) return column;
  }
  return "";
}

function parseElnIds(value: string): string[] {
  return value.trim().split(/\s+/).filter(Boolean);
}

function createBatchAllocator(ids: string[], lastCounts: number[]): () => string {
  let index = 0;
  let batch = lastCounts[0] ?? 0;
  return () => {
    while (index < ids.length && batch >= 999) {
      index += 1;
      batch = lastCounts[index] ?? 0;
    }
    if (index >= ids.length) throw new Error("The selected ELN IDs do not have enough free batch pages.");
    batch += 1;
    return `${ids[index]}${String(batch).padStart(3, "0")}`;
  };
}

function mappedValue(
  fileData: ParsedBulkFile,
  rowIndex: number,
  column: string
): string {
  if (!column || fileData.kind === "nostruct") return "";
  if (fileData.kind === "sdf") return fileData.rows[rowIndex]?.fields[column] ?? "";
  return fileData.rows[rowIndex]?.[column] ?? "";
}

export default function BulkImportScreen({ session, onSessionExpired }: BulkImportScreenProps) {
  const [fileData, setFileData] = useState<ParsedBulkFile | null>(null);
  const [fileName, setFileName] = useState("");
  const [metadata, setMetadata] = useState<UploadMetadata>({ ...emptyMetadata });
  const [chemists, setChemists] = useState<string[]>([]);
  const [compoundTypes, setCompoundTypes] = useState<string[]>([]);
  const [projects, setProjects] = useState<string[]>([]);
  const [suppliers, setSuppliers] = useState<string[]>([]);
  const [solvents, setSolvents] = useState<string[]>([]);
  const [products, setProducts] = useState<string[]>([]);
  const [libraries, setLibraries] = useState<string[]>([]);
  const [elnIdsText, setElnIdsText] = useState("");
  const [externalIdColumn, setExternalIdColumn] = useState("");
  const [batchColumn, setBatchColumn] = useState("");
  const [purityColumn, setPurityColumn] = useState("");
  const [lookupLoading, setLookupLoading] = useState(true);
  const [parsing, setParsing] = useState(false);
  const [uploading, setUploading] = useState(false);
  const [progress, setProgress] = useState({ current: 0, total: 0 });
  const [error, setError] = useState("");
  const [report, setReport] = useState<ImportReport | null>(null);

  useEffect(() => {
    let active = true;
    Promise.all([
      getColumnValues("chemist", session),
      getColumnValues("compound_type", session),
      getColumnValues("project", session),
      getColumnValues("supplier", session),
      getColumnValues("solvent", session),
      getColumnValues("product", session),
      getColumnValues("library_id", session)
    ])
      .then(([loadedChemists, loadedCompoundTypes, loadedProjects, loadedSuppliers, loadedSolvents, loadedProducts, loadedLibraries]) => {
        if (!active) return;
        setChemists(loadedChemists);
        setCompoundTypes(loadedCompoundTypes);
        setProjects(loadedProjects);
        setSuppliers(loadedSuppliers);
        setSolvents(loadedSolvents);
        setProducts(loadedProducts);
        setLibraries(loadedLibraries);
      })
      .catch((loadError: unknown) => {
        if (!active) return;
        setError(errorText(loadError));
        if (loadError instanceof ApiError && loadError.status === 401) onSessionExpired();
      })
      .finally(() => {
        if (active) setLookupLoading(false);
      });
    return () => {
      active = false;
    };
  }, [onSessionExpired, session]);

  async function selectFile(event: ChangeEvent<HTMLInputElement>) {
    const file = event.currentTarget.files?.[0];
    event.currentTarget.value = "";
    if (!file) return;
    setParsing(true);
    setError("");
    setReport(null);
    setFileData(null);
    setFileName(file.name);
    try {
      const parsed = await parseBulkFile(file);
      setFileData(parsed);
      if (parsed.kind !== "nostruct") {
        setExternalIdColumn(chooseDefaultColumn(parsed.columns, ["external_id", "external id", "id"]));
        setBatchColumn(chooseDefaultColumn(parsed.columns, ["supplier_batch", "external batch", "batch"]));
        setPurityColumn(chooseDefaultColumn(parsed.columns, ["purity"]));
      }
    } catch (parseError) {
      setFileName("");
      setError(errorText(parseError));
    } finally {
      setParsing(false);
    }
  }

  function setMetadataValue<Key extends keyof UploadMetadata>(key: Key, value: UploadMetadata[Key]) {
    setMetadata((current) => ({ ...current, [key]: value }));
  }

  function handleError(requestError: unknown) {
    setError(errorText(requestError));
    if (requestError instanceof ApiError && requestError.status === 401) onSessionExpired();
  }

  async function startImport() {
    setError("");
    setReport(null);
    if (!fileData) {
      setError("Choose an SDF or Excel workbook first.");
      return;
    }
    if (lookupLoading) {
      setError("Wait for the registry options to finish loading.");
      return;
    }
    if (!metadata.chemist || !metadata.compoundType || !metadata.project || !metadata.supplier ||
        !metadata.solvent || !metadata.product || !metadata.library || !metadata.ipRights) {
      setError("Complete all registration metadata fields before importing.");
      return;
    }

    const isCtrl = fileData.kind === "nostruct" && metadata.compoundType === "CTRL Compound";
    if (metadata.compoundType === "CTRL Compound" && fileData.kind !== "nostruct") {
      setError("CTRL Compound uploads require the no-structure workbook format.");
      return;
    }

    const elnIds = isCtrl ? [] : parseElnIds(elnIdsText);
    if (!isCtrl && (!elnIds.length || elnIds.some((id) => !/^[a-zA-Z0-9]{6}$/.test(id)))) {
      setError("Enter one or more six-character ELN IDs separated by spaces.");
      return;
    }
    if (new Set(elnIds).size !== elnIds.length) {
      setError("ELN IDs must be unique.");
      return;
    }

    setUploading(true);
    setProgress({ current: 0, total: fileData.rows.length });
    const issues: ImportIssue[] = [];
    let submitted = 0;
    let registered = 0;
    let newCompounds = 0;

    try {
      const lastCounts = isCtrl
        ? []
        : await Promise.all(elnIds.map((id) => getLastBatchOfEln(id, session)));
      const freePages = lastCounts.reduce((available, count) => available + Math.max(0, 999 - count), 0);
      if (!isCtrl && freePages < fileData.rows.length) {
        throw new Error(`Selected ELN IDs have ${freePages} free pages for ${fileData.rows.length} rows.`);
      }

      if (fileData.kind === "nostruct") {
        if (isCtrl) {
          let controlNumber = await getLastCtrlCompound(session);
          for (let index = 0; index < fileData.rows.length; index += 1) {
            const row = fileData.rows[index];
            if (!row) continue;
            controlNumber += 1;
            const result = await postCtrlMol({
              ...uploadMetadataFields(metadata),
              external_id: row.externalId,
              supplier_batch: row.supplierBatch,
              mw: String(row.mw),
              restriction_comment: row.restrictionComment,
              purity: "-1",
              jpage: `CTRL${String(controlNumber).padStart(5, "0")}`,
              iCtrlComp: String(controlNumber)
            }, session);
            if (result.status === 401) throw new ApiError(result.body || "Your session has expired.", 401);
            if (result.ok) {
              submitted += 1;
              registered += 1;
            } else {
              issues.push({ row: String(index + 2), identifier: row.externalId, kind: "error", message: result.body || `HTTP ${result.status}` });
            }
            setProgress({ current: index + 1, total: fileData.rows.length });
          }
        } else {
          const nextPage = createBatchAllocator(elnIds, lastCounts);
          for (let index = 0; index < fileData.rows.length; index += 1) {
            const row = fileData.rows[index];
            if (!row) continue;
            const result = await postNostructMol({
              ...uploadMetadataFields(metadata),
              external_id: row.externalId,
              supplier_batch: row.supplierBatch,
              mw: String(row.mw),
              restriction_comment: row.restrictionComment,
              purity: "-1",
              jpage: nextPage()
            }, session);
            if (result.status === 401) throw new ApiError(result.body || "Your session has expired.", 401);
            if (result.ok) {
              submitted += 1;
              registered += 1;
            } else {
              issues.push({ row: String(index + 2), identifier: row.externalId, kind: "error", message: result.body || `HTTP ${result.status}` });
            }
            setProgress({ current: index + 1, total: fileData.rows.length });
          }
        }
      } else {
        const sequence = await getNextSdfSequence(session);
        const nextPage = createBatchAllocator(elnIds, lastCounts);
        for (let index = 0; index < fileData.rows.length; index += 1) {
          const row = fileData.rows[index];
          if (!row) continue;
          const structure = fileData.kind === "sdf"
            ? fileData.rows[index]?.molfile ?? ""
            : fileData.rows[index]?.Smile ?? "";
          const externalId = mappedValue(fileData, index, externalIdColumn);
          const supplierBatch = mappedValue(fileData, index, batchColumn);
          const purityValue = mappedValue(fileData, index, purityColumn);
          if (!structure.trim()) {
            issues.push({ row: String(index + 1), identifier: externalId, kind: "error", message: "Structure is empty." });
            setProgress({ current: index + 1, total: fileData.rows.length });
            continue;
          }

          const body = {
            ...uploadMetadataFields(metadata),
            [fileData.kind === "sdf" ? "molfile" : "smiles"]: structure,
            external_id: externalId,
            supplier_batch: supplierBatch,
            purity: purityValue || "-1",
            jpage: nextPage(),
            sdfile_sequence: sequence
          };
          const result = await postChemRegMol(body, session);
          if (result.status === 401) throw new ApiError(result.body || "Your session has expired.", 401);
          if (result.ok) {
            submitted += 1;
            const status = result.body.split(";warning=", 1)[0];
            if (status === "newMolecule") newCompounds += 1;
            if (result.body.includes(";warning=")) {
              const warning = result.body.split(";warning=", 2)[1] || "Structure warning";
              if (!warning.includes("Tautomer canonicalization changed connectivity")) {
                issues.push({
                  row: String(index + 1),
                  identifier: externalId,
                  kind: "warning",
                  message: warning
                });
              }
            }
          } else {
            const error = result.body || `HTTP ${result.status}`;
            const message = fileData.kind === "smiles" ? `${error}\nSMILES: ${structure}` : error;
            issues.push({ row: String(index + 1), identifier: externalId, kind: "error", message });
          }
          setProgress({ current: index + 1, total: fileData.rows.length });
        }

        const regnos = await getRegnosFromSdfSequence(sequence, session);
        if (regnos.length !== submitted) {
          issues.push({
            row: "-",
            identifier: `Sequence ${sequence}`,
            kind: "error",
            message: `Expected ${submitted} registered rows, found ${regnos.length}.`
          });
        }
        for (let index = 0; index < regnos.length; index += 1) {
          const regno = regnos[index];
          if (!regno) continue;
          try {
            await registerCompound(regno, session);
            registered += 1;
          } catch (registerError) {
            if (registerError instanceof ApiError && registerError.status === 401) throw registerError;
            issues.push({ row: "-", identifier: regno, kind: "error", message: errorText(registerError) });
          }
          setProgress({ current: fileData.rows.length + index + 1, total: fileData.rows.length + regnos.length });
        }
      }

      setReport({ total: fileData.rows.length, submitted, registered, newCompounds, issues });
    } catch (importError) {
      handleError(importError);
      if (importError instanceof ApiError && importError.status === 401) onSessionExpired();
      if (submitted || issues.length) {
        setReport({ total: fileData.rows.length, submitted, registered, newCompounds, issues });
      }
    } finally {
      setUploading(false);
    }
  }

  function downloadReport() {
    if (!report) return;
    const lines = ["row\tidentifier\ttype\tmessage"];
    for (const issue of report.issues) {
      lines.push([issue.row, issue.identifier, issue.kind, issue.message.replace(/[\r\n\t]+/g, " ")].join("\t"));
    }
    const blob = new Blob([lines.join("\n")], { type: "text/tab-separated-values" });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = "chemreg-import-report.tsv";
    link.click();
    URL.revokeObjectURL(url);
  }

  const rowCount = fileData?.rows.length ?? 0;
  const isCtrlWorkbook = fileData?.kind === "nostruct" && metadata.compoundType === "CTRL Compound";
  const requiredElNCount = fileData && rowCount ? Math.ceil(rowCount / 999) : 0;

  return (
    <section className="bulk-import-screen">
      <header className="bulk-import-heading">
        <div>
          <span className="eyebrow">BULK REGISTRATION / {session.database.toUpperCase()}</span>
          <h2>Import compounds</h2>
        </div>
        <span className="bulk-record-count">{rowCount ? `${rowCount} rows` : "SDF / XLSX"}</span>
      </header>

      {error && <div className="screen-alert" role="alert"><AlertCircle size={17} /><span>{error}</span></div>}

      <div className="bulk-layout">
        <div className="bulk-main">
          <section className="bulk-section">
            <header className="registration-section-heading"><span className="eyebrow">SOURCE FILE</span><span className="section-index">01</span></header>
            <label className="bulk-file-picker">
              <input accept=".sdf,.xlsx,.xls" disabled={parsing || uploading} onChange={(event) => void selectFile(event)} type="file" />
              {parsing ? <LoaderCircle className="spin" size={19} /> : <FileUp size={19} />}
              <span>{parsing ? "Reading file" : fileName || "Choose SDF or XLSX file"}</span>
              <small>{fileData?.kind === "sdf" ? "Structure data file" : fileData?.kind === "smiles" ? "SMILES spreadsheet" : fileData?.kind === "nostruct" ? "No-structure spreadsheet" : "SDF, .xlsx and .xls supported"}</small>
            </label>
            {fileData && fileData.kind !== "nostruct" && (
              <div className="mapping-grid">
                {fileData.kind === "smiles" && <div className="mapping-fixed"><span>STRUCTURE COLUMN</span><strong>Smile</strong></div>}
                <ColumnSelect label="External ID" columns={fileData.columns} value={externalIdColumn} onChange={setExternalIdColumn} />
                <ColumnSelect label="External batch" columns={fileData.columns} value={batchColumn} onChange={setBatchColumn} />
                <ColumnSelect label="Purity" columns={fileData.columns} value={purityColumn} onChange={setPurityColumn} />
              </div>
            )}
          </section>

          <section className="bulk-section">
            <header className="registration-section-heading"><span className="eyebrow">REGISTRATION METADATA</span><span className="section-index">02</span></header>
            <div className="bulk-metadata-grid">
              <SelectField label="Submitter" value={metadata.chemist} options={chemists} onChange={(value) => setMetadataValue("chemist", value)} />
              <SelectField label="Compound type" value={metadata.compoundType} options={compoundTypes} onChange={(value) => {
                setMetadataValue("compoundType", value);
                setElnIdsText(value === "CTRL Compound" && fileData?.kind === "nostruct" ? "CTRL" : "");
              }} />
              <SelectField label="Project" value={metadata.project} options={projects} onChange={(value) => setMetadataValue("project", value)} />
              <SelectField label="Supplier" value={metadata.supplier} options={suppliers} onChange={(value) => setMetadataValue("supplier", value)} />
              <SelectField label="Analytical solvent" value={metadata.solvent} options={solvents} onChange={(value) => setMetadataValue("solvent", value)} />
              <SelectField label="Product" value={metadata.product} options={products} onChange={(value) => setMetadataValue("product", value)} />
              <SelectField label="Library ID" value={metadata.library} options={libraries} onChange={(value) => setMetadataValue("library", value)} />
              <SelectField label="IP rights" value={metadata.ipRights} options={ipRightsOptions} onChange={(value) => setMetadataValue("ipRights", value)} />
            </div>
          </section>

          {fileData && !isCtrlWorkbook && (
            <section className="bulk-section">
              <header className="registration-section-heading"><span className="eyebrow">ELN ALLOCATION</span><span className="section-index">03</span></header>
              <label className="form-field eln-field">
                <span>ELN IDs</span>
                <textarea disabled={metadata.compoundType === "CTRL Compound"} onChange={(event) => setElnIdsText(event.target.value)} placeholder="ABC123 DEF456" rows={2} value={elnIdsText} />
                <small className="field-hint">{requiredElNCount} minimum for {rowCount} rows, up to 999 compounds per ELN ID</small>
              </label>
            </section>
          )}

          <div className="bulk-actions">
            <button className="primary-button bulk-submit" disabled={uploading || parsing || lookupLoading || !fileData} onClick={() => void startImport()} type="button">
              {uploading ? <LoaderCircle className="spin" size={17} /> : <Upload size={17} />}
              <span>{uploading ? "Importing records" : "Start import"}</span>
            </button>
            {report && <span className="bulk-finished"><Check size={15} />Import complete</span>}
          </div>
        </div>

        <aside className="bulk-preview">
          <header className="registration-section-heading"><span className="eyebrow">IMPORT PREVIEW</span><span className="section-index">04</span></header>
          {!fileData ? (
            <div className="bulk-empty"><FileUp size={23} /><span>No source selected</span></div>
          ) : (
            <>
              <div className="bulk-preview-meta"><strong>{fileName}</strong><span>{rowCount} rows / {fileData.kind.toUpperCase()}</span></div>
              <div className="bulk-preview-list">
                {fileData.rows.slice(0, 6).map((_, index) => {
                  let identifier = `Row ${index + 1}`;
                  let structure = "";
                  if (fileData.kind === "nostruct") {
                    const row = fileData.rows[index];
                    identifier = row?.externalId || identifier;
                    structure = `${row?.mw ?? 0} MW`;
                  } else if (fileData.kind === "sdf") {
                    identifier = mappedValue(fileData, index, externalIdColumn) || identifier;
                    structure = "Molfile record";
                  } else {
                    identifier = mappedValue(fileData, index, externalIdColumn) || identifier;
                    structure = fileData.rows[index]?.Smile ?? "";
                  }
                  return (
                    <div className="bulk-preview-row" key={index}>
                      <span className="bulk-row-number">{String(index + 1).padStart(2, "0")}</span>
                      <span className="bulk-row-id">{identifier || `Row ${index + 1}`}</span>
                      <span className="bulk-row-structure">{structure || "Empty structure"}</span>
                    </div>
                  );
                })}
                {rowCount > 6 && <span className="bulk-more-rows">+ {rowCount - 6} more rows</span>}
              </div>
            </>
          )}
          {uploading && (
            <div className="bulk-progress">
              <div><span>Processing</span><span>{progress.current} / {progress.total}</span></div>
              <progress max={Math.max(progress.total, 1)} value={progress.current} />
            </div>
          )}
          <div className="bulk-caveat"><TriangleAlert size={16} /><span>Rows are submitted individually. A failed import can leave earlier rows registered.</span></div>
        </aside>
      </div>

      {report && (
        <section className="import-report">
          <header className="import-report-heading">
            <div><span className="eyebrow">RESULTS</span><h3>Import report</h3></div>
            <button className="secondary-button" onClick={downloadReport} type="button"><Download size={15} /><span>Download report</span></button>
          </header>
          <div className="import-report-stats">
            <ReportStat label="Rows" value={String(report.total)} />
            <ReportStat label="Submitted" value={String(report.submitted)} />
            <ReportStat label="Registered" value={String(report.registered)} />
            <ReportStat label="New compounds" value={String(report.newCompounds)} />
            <ReportStat label="Issues" value={String(report.issues.filter((issue) => issue.kind === "error").length)} />
          </div>
          {report.issues.length > 0 && (
            <div className="import-issue-list">
              {report.issues.map((issue, index) => (
                <div className={issue.kind === "warning" ? "import-issue import-warning" : "import-issue"} key={`${issue.row}-${issue.identifier}-${index}`}>
                  <span>{issue.kind === "warning" ? "WARNING" : `ROW ${issue.row}`}</span>
                  <strong>{issue.identifier || "Import"}</strong>
                  <p>{issue.message}</p>
                </div>
              ))}
            </div>
          )}
        </section>
      )}
    </section>
  );
}

function uploadMetadataFields(metadata: UploadMetadata): Record<string, string> {
  return {
    chemist: metadata.chemist,
    compound_type: metadata.compoundType,
    project: metadata.project,
    source: metadata.supplier,
    solvent: metadata.solvent,
    product: metadata.product,
    library_id: metadata.library,
    ip_rights: metadata.ipRights
  };
}

function SelectField({
  label,
  value,
  options,
  onChange
}: {
  label: string;
  value: string;
  options: string[];
  onChange: (value: string) => void;
}) {
  return (
    <label className="form-field">
      <span>{label}</span>
      <select onChange={(event) => onChange(event.target.value)} value={value}>
        <option value="">Select {label.toLowerCase()}</option>
        {options.map((option) => <option key={option} value={option}>{option}</option>)}
      </select>
    </label>
  );
}

function ColumnSelect({
  label,
  columns,
  value,
  onChange
}: {
  label: string;
  columns: string[];
  value: string;
  onChange: (value: string) => void;
}) {
  return (
    <label className="form-field">
      <span>{label}</span>
      <select onChange={(event) => onChange(event.target.value)} value={value}>
        <option value="">Do not map</option>
        {columns.map((column) => <option key={column} value={column}>{column}</option>)}
      </select>
    </label>
  );
}

function ReportStat({ label, value }: { label: string; value: string }) {
  return <div><span>{label}</span><strong>{value}</strong></div>;
}