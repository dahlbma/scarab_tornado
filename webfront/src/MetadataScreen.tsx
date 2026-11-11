import { useEffect, useState, type FormEvent } from "react";
import {
  AlertCircle,
  Beaker,
  Check,
  FlaskConical,
  Layers3,
  LoaderCircle,
  Plus,
  Truck
} from "lucide-react";
import {
  ApiError,
  createLibrary,
  createSalt,
  createSupplier,
  getCanonicalSmiles,
  getColumnValues,
  type Session
} from "./api";

type MetadataTab = "library" | "supplier" | "salt";

interface MetadataScreenProps {
  session: Session;
  onSessionExpired: () => void;
}

export default function MetadataScreen({ session, onSessionExpired }: MetadataScreenProps) {
  const [activeTab, setActiveTab] = useState<MetadataTab>("library");
  const [suppliers, setSuppliers] = useState<string[]>([]);
  const [libraries, setLibraries] = useState<string[]>([]);
  const [libraryName, setLibraryName] = useState("");
  const [librarySupplier, setLibrarySupplier] = useState("");
  const [restricted, setRestricted] = useState(false);
  const [supplierName, setSupplierName] = useState("");
  const [saltInput, setSaltInput] = useState("");
  const [pendingSalt, setPendingSalt] = useState("");
  const [loading, setLoading] = useState(true);
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState("");
  const [notice, setNotice] = useState("");

  useEffect(() => {
    let active = true;
    Promise.all([
      getColumnValues("supplier", session),
      getColumnValues("library_description", session)
    ])
      .then(([supplierValues, libraryValues]) => {
        if (!active) return;
        setSuppliers(supplierValues);
        setLibraries(libraryValues);
      })
      .catch((loadError: unknown) => {
        if (!active) return;
        setError(errorText(loadError));
        if (loadError instanceof ApiError && loadError.status === 401) onSessionExpired();
      })
      .finally(() => {
        if (active) setLoading(false);
      });
    return () => {
      active = false;
    };
  }, [onSessionExpired, session]);

  function showError(requestError: unknown) {
    setError(errorText(requestError));
    if (requestError instanceof ApiError && requestError.status === 401) onSessionExpired();
  }

  async function saveLibrary(event: FormEvent<HTMLFormElement>) {
    event.preventDefault();
    const name = libraryName.trim();
    if (!name) {
      setError("Enter a library description.");
      return;
    }
    if (libraries.some((value) => value.toLocaleLowerCase() === name.toLocaleLowerCase())) {
      setError(`A library named "${name}" already exists.`);
      return;
    }

    setSaving(true);
    setError("");
    setNotice("");
    try {
      await createLibrary(name, librarySupplier, restricted, session);
      setLibraries((values) => [...values, name].sort((a, b) => a.localeCompare(b)));
      setLibraryName("");
      setLibrarySupplier("");
      setRestricted(false);
      setNotice("Library added.");
    } catch (saveError) {
      showError(saveError);
    } finally {
      setSaving(false);
    }
  }

  async function saveSupplier(event: FormEvent<HTMLFormElement>) {
    event.preventDefault();
    const name = supplierName.trim();
    if (!name) {
      setError("Enter a supplier name.");
      return;
    }
    if (suppliers.some((value) => value.toLocaleLowerCase() === name.toLocaleLowerCase())) {
      setError(`A supplier named "${name}" already exists.`);
      return;
    }

    setSaving(true);
    setError("");
    setNotice("");
    try {
      await createSupplier(name, session);
      setSuppliers((values) => [...values, name].sort((a, b) => a.localeCompare(b)));
      setSupplierName("");
      setNotice("Supplier added.");
    } catch (saveError) {
      showError(saveError);
    } finally {
      setSaving(false);
    }
  }

  async function checkSalt(event: FormEvent<HTMLFormElement>) {
    event.preventDefault();
    const smiles = saltInput.trim();
    if (!smiles) {
      setError("Enter a SMILES string.");
      return;
    }

    setSaving(true);
    setError("");
    setNotice("");
    setPendingSalt("");
    try {
      const canonical = await getCanonicalSmiles(smiles, session);
      if (canonical.exists) {
        setError(`This salt is already registered as ${canonical.smiles}.`);
      } else if (!canonical.smiles) {
        setError("The server could not canonicalize this SMILES string.");
      } else {
        setPendingSalt(canonical.smiles);
      }
    } catch (checkError) {
      showError(checkError);
    } finally {
      setSaving(false);
    }
  }

  async function confirmSalt() {
    if (!pendingSalt) return;
    setSaving(true);
    setError("");
    try {
      await createSalt(pendingSalt, session);
      setNotice("Salt added.");
      setSaltInput("");
      setPendingSalt("");
    } catch (saveError) {
      showError(saveError);
    } finally {
      setSaving(false);
    }
  }

  return (
    <section className="metadata-screen">
      <div className="metadata-screen-grid">
        <aside className="metadata-nav-panel">
          <span className="eyebrow">REGISTRY SETTINGS</span>
          <h2>Metadata</h2>
          <div className="metadata-tabs" role="tablist" aria-label="Metadata type">
            <MetadataTabButton active={activeTab === "library"} icon={<Layers3 size={17} />} label="Libraries" onClick={() => setActiveTab("library")} />
            <MetadataTabButton active={activeTab === "supplier"} icon={<Truck size={17} />} label="Suppliers" onClick={() => setActiveTab("supplier")} />
            <MetadataTabButton active={activeTab === "salt"} icon={<Beaker size={17} />} label="Salts" onClick={() => setActiveTab("salt")} />
          </div>
          <div className="metadata-nav-foot">
            <FlaskConical size={16} />
            <span>{session.database.toUpperCase()} DATABASE</span>
          </div>
        </aside>

        <div className="metadata-main-panel">
          <header className="metadata-main-heading">
            <div>
              <span className="eyebrow">CHEMREG / METADATA</span>
              <h2>{activeTab === "library" ? "Add library" : activeTab === "supplier" ? "Add supplier" : "Add salt"}</h2>
            </div>
            <span className="metadata-total">
              {activeTab === "library" ? `${libraries.length} libraries` : activeTab === "supplier" ? `${suppliers.length} suppliers` : "SMILES"}
            </span>
          </header>

          {error && <div className="screen-alert" role="alert"><AlertCircle size={17} /><span>{error}</span></div>}
          {notice && <div className="screen-notice" role="status"><Check size={16} /><span>{notice}</span></div>}

          {loading ? (
            <div className="screen-loading"><LoaderCircle className="spin" size={20} />Loading metadata</div>
          ) : activeTab === "library" ? (
            <form className="metadata-form" onSubmit={saveLibrary}>
              <label className="form-field">
                <span>Library description</span>
                <input autoComplete="off" onChange={(event) => setLibraryName(event.target.value)} placeholder="Enter a library name" value={libraryName} />
              </label>
              <label className="form-field">
                <span>Supplier</span>
                <select onChange={(event) => setLibrarySupplier(event.target.value)} value={librarySupplier}>
                  <option value="">Select supplier</option>
                  {suppliers.map((supplier) => <option key={supplier} value={supplier}>{supplier}</option>)}
                </select>
              </label>
              <label className="metadata-toggle">
                <input checked={restricted} onChange={(event) => setRestricted(event.target.checked)} type="checkbox" />
                <span className="toggle-track" aria-hidden="true"><i /></span>
                <span>Restricted library</span>
              </label>
              <button className="primary-button metadata-submit" disabled={saving} type="submit">
                {saving ? <LoaderCircle className="spin" size={17} /> : <Plus size={17} />}
                <span>{saving ? "Adding library" : "Add library"}</span>
              </button>
              <div className="metadata-existing">
                <span className="eyebrow">EXISTING LIBRARIES</span>
                <div className="metadata-chip-list">
                  {libraries.length ? libraries.map((library) => <span className="metadata-chip" key={library}>{library}</span>) : <span className="metadata-muted">No libraries returned.</span>}
                </div>
              </div>
            </form>
          ) : activeTab === "supplier" ? (
            <form className="metadata-form" onSubmit={saveSupplier}>
              <label className="form-field">
                <span>Supplier name</span>
                <input autoComplete="off" onChange={(event) => setSupplierName(event.target.value)} placeholder="Enter a supplier name" value={supplierName} />
              </label>
              <button className="primary-button metadata-submit" disabled={saving} type="submit">
                {saving ? <LoaderCircle className="spin" size={17} /> : <Plus size={17} />}
                <span>{saving ? "Adding supplier" : "Add supplier"}</span>
              </button>
              <div className="metadata-existing">
                <span className="eyebrow">EXISTING SUPPLIERS</span>
                <div className="metadata-chip-list">
                  {suppliers.length ? suppliers.map((supplier) => <span className="metadata-chip" key={supplier}>{supplier}</span>) : <span className="metadata-muted">No suppliers returned.</span>}
                </div>
              </div>
            </form>
          ) : (
            <div className="metadata-form">
              <form className="salt-form" onSubmit={checkSalt}>
                <label className="form-field">
                  <span>SMILES</span>
                  <textarea autoComplete="off" onChange={(event) => { setSaltInput(event.target.value); setPendingSalt(""); }} placeholder="Enter salt SMILES" rows={4} value={saltInput} />
                </label>
                <button className="primary-button metadata-submit" disabled={saving} type="submit">
                  {saving ? <LoaderCircle className="spin" size={17} /> : <Beaker size={17} />}
                  <span>{saving ? "Checking structure" : "Check SMILES"}</span>
                </button>
              </form>
              {pendingSalt && (
                <div className="salt-confirmation">
                  <span className="eyebrow">CANONICAL SMILES</span>
                  <code>{pendingSalt}</code>
                  <p>This structure is not in the salt registry. Add it?</p>
                  <div className="salt-confirmation-actions">
                    <button className="secondary-button" disabled={saving} onClick={() => setPendingSalt("")} type="button">Cancel</button>
                    <button className="primary-button" disabled={saving} onClick={() => void confirmSalt()} type="button">
                      {saving ? <LoaderCircle className="spin" size={16} /> : <Plus size={16} />}
                      <span>Add salt</span>
                    </button>
                  </div>
                </div>
              )}
            </div>
          )}
        </div>
      </div>
    </section>
  );
}

function MetadataTabButton({
  active,
  icon,
  label,
  onClick
}: {
  active: boolean;
  icon: React.ReactNode;
  label: string;
  onClick: () => void;
}) {
  return (
    <button aria-selected={active} className={active ? "metadata-tab metadata-tab-active" : "metadata-tab"} onClick={onClick} role="tab" type="button">
      {icon}<span>{label}</span>
    </button>
  );
}

function errorText(error: unknown): string {
  return error instanceof Error ? error.message : "The request could not be completed.";
}