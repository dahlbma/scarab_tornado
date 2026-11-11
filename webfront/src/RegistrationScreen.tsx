import { useEffect, useRef, useState, type ChangeEvent } from "react";
import {
  AlertCircle,
  ArrowDownToLine,
  Check,
  FileUp,
  FlaskConical,
  LoaderCircle,
  Save,
  X
} from "lucide-react";
import {
  ApiError,
  createDraftRegno,
  createMolImageUrl,
  deleteDraftRegno,
  getColumnValues,
  getMolfile,
  getRegistrationRecord,
  registerCompound,
  updateCompoundValue,
  updateRegistrationBatch,
  uploadMolfile,
  type RegistrationRecord,
  type Session
} from "./api";

interface RegistrationForm {
  externalId: string;
  externalBatch: string;
  batch: string;
  chemist: string;
  project: string;
  compoundType: string;
  product: string;
  library: string;
  solvent: string;
  ipRights: string;
  purity: string;
  chromText: string;
  chromPurity: string;
  nmrText: string;
  nmrPurity: string;
  msText: string;
  msPurity: string;
  comments: string;
}

interface RegistrationScreenProps {
  requestedRegno: string | null;
  session: Session;
  onCancel: () => void;
  onComplete: (regno: string) => void;
  onSessionExpired: () => void;
}

const emptyForm: RegistrationForm = {
  externalId: "",
  externalBatch: "",
  batch: "",
  chemist: "",
  project: "",
  compoundType: "",
  product: "",
  library: "",
  solvent: "",
  ipRights: "",
  purity: "",
  chromText: "",
  chromPurity: "",
  nmrText: "",
  nmrPurity: "",
  msText: "",
  msPurity: "",
  comments: ""
};

const updateColumns: Array<[keyof RegistrationForm, string]> = [
  ["externalId", "EXTERNAL_ID"],
  ["externalBatch", "SUPPLIER_BATCH"],
  ["chemist", "CHEMIST"],
  ["project", "project"],
  ["compoundType", "compound_type"],
  ["product", "product"],
  ["library", "library_id"],
  ["solvent", "solvent"],
  ["ipRights", "ip_rights"],
  ["purity", "PURITY"],
  ["chromText", "CHROM_TEXT"],
  ["chromPurity", "CHROM_PURITY"],
  ["nmrText", "NMR_TEXT"],
  ["nmrPurity", "NMR_PURITY"],
  ["msText", "MS_TEXT"],
  ["msPurity", "MS_PURITY"],
  ["comments", "COMMENTS"]
];

const ipRightsOptions = ["External rights", "LCBKI", "Commercial"];
const batchPattern = /^[a-zA-Z]{2}[0-9]{7}$/;

function formFromRecord(record: RegistrationRecord): RegistrationForm {
  return {
    externalId: record.externalId,
    externalBatch: record.externalBatch,
    batch: record.jpage,
    chemist: record.chemist,
    project: record.project,
    compoundType: record.compoundType,
    product: record.product,
    library: record.library,
    solvent: record.solvent,
    ipRights: record.ipRights,
    purity: record.purity,
    chromText: record.chromText,
    chromPurity: record.chromPurity,
    nmrText: record.nmrText,
    nmrPurity: record.nmrPurity,
    msText: record.msText,
    msPurity: record.msPurity,
    comments: record.comments
  };
}

function errorText(error: unknown): string {
  return error instanceof Error ? error.message : "The request could not be completed.";
}

export default function RegistrationScreen({
  requestedRegno,
  session,
  onCancel,
  onComplete,
  onSessionExpired
}: RegistrationScreenProps) {
  const [regno, setRegno] = useState(requestedRegno ?? "");
  const [record, setRecord] = useState<RegistrationRecord | null>(null);
  const [form, setForm] = useState<RegistrationForm>({ ...emptyForm });
  const [chemists, setChemists] = useState<string[]>([]);
  const [projects, setProjects] = useState<string[]>([]);
  const [compoundTypes, setCompoundTypes] = useState<string[]>([]);
  const [products, setProducts] = useState<string[]>([]);
  const [libraries, setLibraries] = useState<string[]>([]);
  const [solvents, setSolvents] = useState<string[]>([]);
  const [hasStructure, setHasStructure] = useState(false);
  const [imageUrl, setImageUrl] = useState("");
  const [structureName, setStructureName] = useState("");
  const [loading, setLoading] = useState(true);
  const [saving, setSaving] = useState(false);
  const [uploading, setUploading] = useState(false);
  const [registering, setRegistering] = useState(false);
  const [createdDraft, setCreatedDraft] = useState(false);
  const [persistedChanges, setPersistedChanges] = useState(false);
  const [dirty, setDirty] = useState(false);
  const [error, setError] = useState("");
  const [notice, setNotice] = useState("");
  const draftRequest = useRef<Promise<string> | null>(null);

  useEffect(() => {
    let active = true;

    async function initialize() {
      setLoading(true);
      setError("");
      try {
        const nextRegno = requestedRegno ?? await (draftRequest.current ??= createDraftRegno(session));
        if (!active) return;
        setRegno(nextRegno);
        if (!requestedRegno) setCreatedDraft(true);

        const [loadedChemists, loadedProjects, loadedCompoundTypes, loadedProducts, loadedLibraries, loadedSolvents] = await Promise.all([
          getColumnValues("chemist", session),
          getColumnValues("project", session),
          getColumnValues("compound_type", session),
          getColumnValues("product", session),
          getColumnValues("library_id", session),
          getColumnValues("solvent", session)
        ]);
        if (!active) return;
        setChemists(loadedChemists);
        setProjects(loadedProjects);
        setCompoundTypes(loadedCompoundTypes);
        setProducts(loadedProducts);
        setLibraries(loadedLibraries);
        setSolvents(loadedSolvents);

        if (requestedRegno) {
          const [currentRecord, molfile, image] = await Promise.all([
            getRegistrationRecord(nextRegno, session),
            getMolfile(nextRegno, session),
            createMolImageUrl(nextRegno, session).catch(() => "")
          ]);
          if (!active) return;
          setRecord(currentRecord);
          setForm(formFromRecord(currentRecord));
          setHasStructure(Boolean(molfile.trim()));
          setImageUrl(image);
        }
      } catch (loadError) {
        if (!active) return;
        setError(errorText(loadError));
        if (loadError instanceof ApiError && loadError.status === 401) onSessionExpired();
      } finally {
        if (active) setLoading(false);
      }
    }

    void initialize();
    return () => {
      active = false;
    };
  }, [onSessionExpired, requestedRegno, session]);

  function setField<Key extends keyof RegistrationForm>(key: Key, value: RegistrationForm[Key]) {
    setForm((current) => ({ ...current, [key]: value }));
    setDirty(true);
    setNotice("");
  }

  function handleError(requestError: unknown) {
    setError(errorText(requestError));
    if (requestError instanceof ApiError && requestError.status === 401) onSessionExpired();
  }

  async function saveDraft(): Promise<boolean> {
    if (!regno) return false;
    setSaving(true);
    setError("");
    setNotice("");
    try {
      for (const [key, column] of updateColumns) {
        const value = form[key] || (key === "purity" ? "-1" : "");
        await updateCompoundValue(regno, column, value, session);
      }

      const validBatch = batchPattern.test(form.batch.trim());
      const currentCompoundId = record?.compoundId.trim().toLowerCase();
      const registered = Boolean(requestedRegno && currentCompoundId && !["none", "null"].includes(currentCompoundId));
      if (validBatch && !registered) {
        await updateRegistrationBatch(regno, form.batch.trim(), session);
      }

      setPersistedChanges(true);
      setDirty(false);
      setNotice(validBatch || registered
        ? "Changes saved."
        : "Other fields saved. Enter a batch such as AB1234567 before saving the batch.");
      return true;
    } catch (saveError) {
      handleError(saveError);
      return false;
    } finally {
      setSaving(false);
    }
  }

  async function chooseMolfile(event: ChangeEvent<HTMLInputElement>) {
    const file = event.currentTarget.files?.[0];
    event.currentTarget.value = "";
    if (!file || !regno) return;
    if (!file.name.toLowerCase().endsWith(".mol")) {
      setError("Choose a molfile with the .mol extension.");
      return;
    }

    setUploading(true);
    setError("");
    setNotice("");
    try {
      await uploadMolfile(regno, file, session);
      const [updatedRecord, updatedImage] = await Promise.all([
        getRegistrationRecord(regno, session),
        createMolImageUrl(regno, session)
      ]);
      setRecord(updatedRecord);
      setImageUrl(updatedImage);
      setHasStructure(true);
      setStructureName(file.name);
      setPersistedChanges(true);
      setNotice("Structure uploaded and molecular properties refreshed.");
    } catch (uploadError) {
      handleError(uploadError);
    } finally {
      setUploading(false);
    }
  }

  async function downloadMolfile() {
    if (!regno) return;
    try {
      const molfile = await getMolfile(regno, session);
      if (!molfile.trim()) throw new Error("No molfile is stored for this registration.");
      const url = URL.createObjectURL(new Blob([molfile], { type: "chemical/x-mdl-molfile" }));
      const link = document.createElement("a");
      link.href = url;
      link.download = `${regno}.mol`;
      link.click();
      URL.revokeObjectURL(url);
    } catch (downloadError) {
      handleError(downloadError);
    }
  }

  async function finishRegistration() {
    if (!canRegister) {
      setError("Choose a submitter, project, compound type, product, valid batch, IP rights, and molfile first.");
      return;
    }
    const saved = await saveDraft();
    if (!saved) return;
    setRegistering(true);
    setError("");
    try {
      await registerCompound(regno, session);
      onComplete(regno);
    } catch (registerError) {
      handleError(registerError);
    } finally {
      setRegistering(false);
    }
  }

  async function cancelRegistration() {
    if (createdDraft && !persistedChanges && regno) {
      try {
        await deleteDraftRegno(regno, session);
      } catch (deleteError) {
        handleError(deleteError);
        return;
      }
    }
    onCancel();
  }

  const batchValid = batchPattern.test(form.batch.trim());
  const currentCompoundId = record?.compoundId.trim().toLowerCase();
  const isRegistered = Boolean(requestedRegno && currentCompoundId && !["none", "null"].includes(currentCompoundId));
  const canRegister = Boolean(
    regno && form.chemist && form.project && form.compoundType && form.product &&
    form.ipRights && batchValid && hasStructure && !isRegistered && !saving &&
    !uploading && !registering
  );

  if (loading) {
    return <div className="screen-loading"><LoaderCircle className="spin" size={20} />Loading registration form</div>;
  }

  return (
    <section className="registration-screen">
      <div className="registration-toolbar">
        <div>
          <span className="eyebrow">{requestedRegno ? "EDIT REGISTRATION" : "NEW REGISTRATION"}</span>
          <h2>Reg no. <span className="mono-value">{regno || "-"}</span></h2>
        </div>
        <div className="registration-toolbar-actions">
          {isRegistered && <span className="registered-badge"><Check size={14} />COMPOUND {record?.compoundId}</span>}
          {requestedRegno && hasStructure && (
            <button className="secondary-button" onClick={() => void downloadMolfile()} type="button">
              <ArrowDownToLine size={16} /><span>Download molfile</span>
            </button>
          )}
        </div>
      </div>

      {error && <div className="screen-alert" role="alert"><AlertCircle size={17} /><span>{error}</span></div>}
      {notice && <div className="screen-notice" role="status"><Check size={16} /><span>{notice}</span></div>}

      <div className="registration-layout">
        <div className="registration-main">
          <section className="registration-section">
            <header className="registration-section-heading">
              <span className="eyebrow">IDENTIFIERS / OWNERSHIP</span>
              <span className="section-index">01</span>
            </header>
            <div className="registration-fields registration-fields-three">
              <SelectField label="Submitter" value={form.chemist} options={chemists} onChange={(value) => setField("chemist", value)} />
              <SelectField label="Project" value={form.project} options={projects} onChange={(value) => setField("project", value)} />
              <SelectField label="Compound type" value={form.compoundType} options={compoundTypes} onChange={(value) => setField("compoundType", value)} />
              <SelectField label="Product" value={form.product} options={products} onChange={(value) => setField("product", value)} />
              <SelectField label="Library ID" value={form.library} options={libraries} onChange={(value) => setField("library", value)} />
              <SelectField label="Analytical solvent" value={form.solvent} options={solvents} onChange={(value) => setField("solvent", value)} />
              <TextField label="External ID" value={form.externalId} onChange={(value) => setField("externalId", value)} />
              <TextField label="External batch" value={form.externalBatch} onChange={(value) => setField("externalBatch", value)} />
              <label className="form-field">
                <span>ELN batch</span>
                <input
                  aria-invalid={Boolean(form.batch) && !batchValid}
                  disabled={isRegistered}
                  maxLength={9}
                  onChange={(event) => setField("batch", event.target.value)}
                  placeholder="AB1234567"
                  value={form.batch}
                />
                <small className={form.batch && !batchValid ? "field-error" : "field-hint"}>Two letters followed by seven digits</small>
              </label>
              <SelectField label="IP rights" value={form.ipRights} options={ipRightsOptions} onChange={(value) => setField("ipRights", value)} />
              <label className="form-field">
                <span>Purity (%)</span>
                <input
                  inputMode="numeric"
                  max="100"
                  min="0"
                  onChange={(event) => setField("purity", event.target.value)}
                  placeholder="-"
                  type="number"
                  value={form.purity}
                />
              </label>
            </div>
          </section>

          <section className="registration-section">
            <header className="registration-section-heading">
              <span className="eyebrow">ANALYTICAL DATA</span>
              <span className="section-index">02</span>
            </header>
            <div className="analysis-grid">
              <AnalysisField label="Chromatography" text={form.chromText} purity={form.chromPurity} onTextChange={(value) => setField("chromText", value)} onPurityChange={(value) => setField("chromPurity", value)} />
              <AnalysisField label="NMR" text={form.nmrText} purity={form.nmrPurity} onTextChange={(value) => setField("nmrText", value)} onPurityChange={(value) => setField("nmrPurity", value)} />
              <AnalysisField label="Mass spectrometry" text={form.msText} purity={form.msPurity} onTextChange={(value) => setField("msText", value)} onPurityChange={(value) => setField("msPurity", value)} />
              <label className="form-field comments-field">
                <span>Comments</span>
                <textarea onChange={(event) => setField("comments", event.target.value)} rows={4} value={form.comments} />
              </label>
            </div>
          </section>
        </div>

        <aside className="registration-structure">
          <div className="registration-section-heading">
            <span className="eyebrow">MOLECULAR STRUCTURE</span>
            <span className="section-index">03</span>
          </div>
          <div className="registration-molview">
            {imageUrl ? <img alt={`Molecular structure for ${regno}`} src={imageUrl} /> : <div className="molview-empty"><FlaskConical size={28} /><span>{hasStructure ? "Structure preview unavailable" : "No molfile loaded"}</span></div>}
          </div>
          <div className="molecule-properties">
            <Property label="Formula" value={record?.formula} />
            <Property label="Average mass" value={record?.averageMass} />
            <Property label="Monoisotopic mass" value={record?.monoisotopicMass} />
            <Property label="CHNS" value={record?.chns} />
          </div>
          <label className={isRegistered ? "file-picker file-picker-disabled" : "file-picker"}>
            <input accept=".mol" disabled={isRegistered || uploading || loading} onChange={(event) => void chooseMolfile(event)} type="file" />
            {uploading ? <LoaderCircle className="spin" size={17} /> : <FileUp size={17} />}
            <span>{uploading ? "Uploading molfile" : hasStructure ? "Replace molfile" : "Choose molfile"}</span>
          </label>
          <span className="file-name">{structureName || (hasStructure ? "Structure stored for this registration" : "Accepted format: .mol")}</span>
          {isRegistered && <span className="field-hint">Structure and ELN batch are locked after compound registration.</span>}
        </aside>
      </div>

      <footer className="registration-actions">
        <button className="secondary-button" disabled={saving || uploading || registering} onClick={() => void cancelRegistration()} type="button">
          <X size={16} /><span>Cancel</span>
        </button>
        <div className="registration-actions-right">
          {dirty && <span className="unsaved-note">Unsaved changes</span>}
          <button className="secondary-button" disabled={saving || uploading || registering || !regno} onClick={() => void saveDraft()} type="button">
            {saving ? <LoaderCircle className="spin" size={16} /> : <Save size={16} />}
            <span>{saving ? "Saving" : "Save draft"}</span>
          </button>
          {!isRegistered && (
            <button className="primary-button" disabled={!canRegister} onClick={() => void finishRegistration()} type="button">
              {registering ? <LoaderCircle className="spin" size={17} /> : <Check size={17} />}
              <span>{registering ? "Registering" : "Register compound"}</span>
            </button>
          )}
        </div>
      </footer>
    </section>
  );
}

function TextField({ label, value, onChange }: { label: string; value: string; onChange: (value: string) => void }) {
  return (
    <label className="form-field">
      <span>{label}</span>
      <input onChange={(event) => onChange(event.target.value)} value={value} />
    </label>
  );
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

function AnalysisField({
  label,
  text,
  purity,
  onTextChange,
  onPurityChange
}: {
  label: string;
  text: string;
  purity: string;
  onTextChange: (value: string) => void;
  onPurityChange: (value: string) => void;
}) {
  return (
    <div className="analysis-field">
      <label className="form-field">
        <span>{label}</span>
        <textarea onChange={(event) => onTextChange(event.target.value)} rows={4} value={text} />
      </label>
      <label className="form-field purity-field">
        <span>Purity (%)</span>
        <input inputMode="numeric" max="100" min="0" onChange={(event) => onPurityChange(event.target.value)} type="number" value={purity} />
      </label>
    </div>
  );
}

function Property({ label, value }: { label: string; value?: string }) {
  return <div><span>{label}</span><strong>{value || "-"}</strong></div>;
}