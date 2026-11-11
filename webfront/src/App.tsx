import { useCallback, useEffect, useState, type FormEvent } from "react";
import {
  Activity,
  Atom,
  ChevronLeft,
  ChevronRight,
  CircleAlert,
  Database,
  FlaskConical,
  KeyRound,
  LoaderCircle,
  LogOut,
  Pencil,
  RotateCcw,
  ScanSearch,
  Search,
  UserRound
} from "lucide-react";
import {
  ApiError,
  createMolImageUrl,
  getColumnValues,
  getCompoundRecord,
  getDatabases,
  login,
  searchCompounds,
  type CompoundRecord,
  type SearchMode,
  type Session
} from "./api";
import RegistrationScreen from "./RegistrationScreen";
import MetadataScreen from "./MetadataScreen";
import BulkImportScreen from "./BulkImportScreen";

const searchModes: Array<{ value: SearchMode; label: string }> = [
  { value: "regno", label: "Reg no." },
  { value: "submitter", label: "Submitter" },
  { value: "compoundId", label: "Compound ID" },
  { value: "batch", label: "Batch" }
];

const searchPlaceholders: Record<SearchMode, string> = {
  regno: "Enter a registration number",
  submitter: "Choose or enter a submitter",
  compoundId: "Enter a compound ID",
  batch: "Enter a batch or ELN prefix"
};

function errorText(error: unknown): string {
  return error instanceof Error ? error.message : "The request could not be completed.";
}

export default function App() {
  const [databases, setDatabases] = useState(["Test", "Live"]);
  const [backendStatus, setBackendStatus] = useState<"checking" | "ready" | "offline">("checking");
  const [session, setSession] = useState<Session | null>(null);
  const [loginError, setLoginError] = useState("");
  const [loginBusy, setLoginBusy] = useState(false);

  useEffect(() => {
    let active = true;
    getDatabases()
      .then((availableDatabases) => {
        if (!active) return;
        if (availableDatabases.length) setDatabases(availableDatabases);
        setBackendStatus("ready");
      })
      .catch(() => {
        if (active) setBackendStatus("offline");
      });
    return () => {
      active = false;
    };
  }, []);

  async function handleLogin(username: string, password: string, database: string) {
    setLoginBusy(true);
    setLoginError("");
    try {
      setSession(await login(username, password, database));
    } catch (error) {
      setLoginError(errorText(error));
    } finally {
      setLoginBusy(false);
    }
  }

  const handleLogout = useCallback((message = "") => {
    setSession(null);
    setLoginError(message);
  }, []);

  if (!session) {
    return (
      <LoginScreen
        databases={databases}
        backendStatus={backendStatus}
        busy={loginBusy}
        error={loginError}
        onLogin={handleLogin}
      />
    );
  }

  return <Workspace session={session} onLogout={handleLogout} />;
}

interface LoginScreenProps {
  databases: string[];
  backendStatus: "checking" | "ready" | "offline";
  busy: boolean;
  error: string;
  onLogin: (username: string, password: string, database: string) => Promise<void>;
}

function LoginScreen({ databases, backendStatus, busy, error, onLogin }: LoginScreenProps) {
  const [username, setUsername] = useState("");
  const [password, setPassword] = useState("");
  const [database, setDatabase] = useState("Test");

  useEffect(() => {
    setDatabase((current) => {
      if (databases.includes(current)) return current;
      if (databases.includes("Test")) return "Test";
      return databases[0] ?? "Test";
    });
  }, [databases]);

  function submit(event: FormEvent<HTMLFormElement>) {
    event.preventDefault();
    void onLogin(username.trim(), password, database);
  }

  return (
    <main className="auth-page">
      <div className="auth-grid" aria-hidden="true" />
      <section className="auth-frame" aria-label="ChemReg sign in">
        <aside className="auth-identity">
          <div className="brand-lockup">
            <span className="brand-mark"><FlaskConical size={21} strokeWidth={1.8} /></span>
            <span>ChemReg</span>
          </div>
          <div className="identity-title">
            <span className="eyebrow">CHEMISTRY / REGISTRATION</span>
            <h1>Structure records<br />&amp; metadata</h1>
            <div className="identity-rule"><i /><i /><i /></div>
          </div>
          <div className="identity-foot">
            <span>SCILIFELAB</span>
            <span className="identity-dot" />
            <span>REGISTRY WORKSPACE</span>
          </div>
        </aside>

        <div className="auth-content">
          <div className="auth-heading">
            <span className="eyebrow">SECURE ACCESS</span>
            <h2>Sign in</h2>
            <p>Use your ChemReg database credentials.</p>
          </div>

          <form className="login-form" onSubmit={submit}>
            <label className="form-field">
              <span>Username</span>
              <span className="input-wrap">
                <UserRound size={17} aria-hidden="true" />
                <input
                  autoComplete="username"
                  autoFocus
                  onChange={(event) => setUsername(event.target.value)}
                  placeholder="Username"
                  required
                  value={username}
                />
              </span>
            </label>

            <label className="form-field">
              <span>Password</span>
              <span className="input-wrap">
                <KeyRound size={17} aria-hidden="true" />
                <input
                  autoComplete="current-password"
                  onChange={(event) => setPassword(event.target.value)}
                  placeholder="Password"
                  required
                  type="password"
                  value={password}
                />
              </span>
            </label>

            <label className="form-field">
              <span>Database</span>
              <span className="input-wrap select-wrap">
                <Database size={17} aria-hidden="true" />
                <select onChange={(event) => setDatabase(event.target.value)} value={database}>
                  {databases.map((name) => <option key={name} value={name}>{name}</option>)}
                </select>
              </span>
            </label>

            {error && (
              <div className="form-alert" role="alert">
                <CircleAlert size={17} />
                <span>{error}</span>
              </div>
            )}

            {backendStatus === "offline" && !error && (
              <div className="form-alert form-alert-muted" role="status">
                <CircleAlert size={17} />
                <span>Cannot reach the ChemReg server at the configured API address.</span>
              </div>
            )}

            <button className="primary-button login-button" disabled={busy} type="submit">
              {busy ? <LoaderCircle className="spin" size={18} /> : <span>Sign in</span>}
              {!busy && <ChevronRight size={18} />}
            </button>
          </form>

          <div className="auth-footer">
            <span className={`status-led status-${backendStatus}`} />
            <span>{backendStatus === "ready" ? "SERVER AVAILABLE" : backendStatus === "checking" ? "CHECKING SERVER" : "SERVER OFFLINE"}</span>
            <span className="auth-footer-version">CHEMREG WEB</span>
          </div>
        </div>
      </section>
    </main>
  );
}

interface WorkspaceProps {
  session: Session;
  onLogout: (message?: string) => void;
}

function Workspace({ session, onLogout }: WorkspaceProps) {
  const [activeView, setActiveView] = useState<"search" | "registration" | "metadata" | "bulk">("search");
  const [requestedRegno, setRequestedRegno] = useState<string | null>(null);
  const [mode, setMode] = useState<SearchMode>("regno");
  const [searchValue, setSearchValue] = useState("");
  const [chemists, setChemists] = useState<string[]>([]);
  const [results, setResults] = useState<string[]>([]);
  const [selectedRegno, setSelectedRegno] = useState<string | null>(null);
  const [searched, setSearched] = useState(false);
  const [searchBusy, setSearchBusy] = useState(false);
  const [recordBusy, setRecordBusy] = useState(false);
  const [searchError, setSearchError] = useState("");
  const [recordError, setRecordError] = useState("");
  const [chemistError, setChemistError] = useState("");
  const [recordRevision, setRecordRevision] = useState(0);
  const [record, setRecord] = useState<CompoundRecord | null>(null);
  const [molImageUrl, setMolImageUrl] = useState("");
  const [molImageFailed, setMolImageFailed] = useState(false);
  const expireSession = useCallback(() => {
    onLogout("Your session has expired. Sign in again.");
  }, [onLogout]);

  useEffect(() => {
    let active = true;
    getColumnValues("chemist", session)
      .then((values) => {
        if (active) setChemists(values);
      })
      .catch((error: unknown) => {
        if (!active) return;
        if (error instanceof ApiError && error.status === 401) {
          onLogout("Your session has expired. Sign in again.");
        } else {
          setChemistError("Submitter suggestions are unavailable.");
        }
      });
    return () => {
      active = false;
    };
  }, [onLogout, session]);

  useEffect(() => {
    if (!selectedRegno) {
      setRecord(null);
      setMolImageUrl("");
      setMolImageFailed(false);
      return;
    }

    let active = true;
    setRecord(null);
    setRecordError("");
    setRecordBusy(true);
    setMolImageUrl("");
    setMolImageFailed(false);

    getCompoundRecord(selectedRegno, session)
      .then((details) => {
        if (active) setRecord(details);
      })
      .catch((error: unknown) => {
        if (!active) return;
        setRecordError(errorText(error));
        if (error instanceof ApiError && error.status === 401) {
          onLogout("Your session has expired. Sign in again.");
        }
      })
      .finally(() => {
        if (active) setRecordBusy(false);
      });

    createMolImageUrl(selectedRegno, session)
      .then((imageUrl) => {
        if (active) setMolImageUrl(imageUrl);
      })
      .catch((error: unknown) => {
        if (active && error instanceof ApiError && error.status === 401) {
          onLogout("Your session has expired. Sign in again.");
        }
      });

    return () => {
      active = false;
    };
  }, [onLogout, recordRevision, selectedRegno, session]);

  async function runSearch(event: FormEvent<HTMLFormElement>) {
    event.preventDefault();
    const value = searchValue.trim();
    if (!value) {
      setSearchError(mode === "submitter" ? "Enter a submitter name." : "Enter a search value.");
      return;
    }

    setSearchBusy(true);
    setSearchError("");
    setRecordError("");
    setResults([]);
    setSelectedRegno(null);
    setSearched(true);
    try {
      const hits = await searchCompounds(mode, value, session);
      setResults(hits);
      setSelectedRegno(hits[0] ?? null);
    } catch (error) {
      setSearchError(errorText(error));
      if (error instanceof ApiError && error.status === 401) {
        onLogout("Your session has expired. Sign in again.");
      }
    } finally {
      setSearchBusy(false);
    }
  }

  function moveRecord(offset: number) {
    if (!selectedRegno) return;
    const currentIndex = results.indexOf(selectedRegno);
    const nextRegno = results[currentIndex + offset];
    if (nextRegno) setSelectedRegno(nextRegno);
  }

  const selectedIndex = selectedRegno ? results.indexOf(selectedRegno) : -1;

  return (
    <main className="workspace">
      <header className="topbar">
        <div className="brand-lockup topbar-brand">
          <span className="brand-mark"><FlaskConical size={19} strokeWidth={1.8} /></span>
          <span>ChemReg</span>
        </div>
        <nav className="workspace-nav" aria-label="Primary navigation">
          <button
            aria-pressed={activeView === "search"}
            className={activeView === "search" ? "workspace-nav-item workspace-nav-item-active" : "workspace-nav-item"}
            onClick={() => setActiveView("search")}
            type="button"
          >Search</button>
          <button
            aria-pressed={activeView === "registration"}
            className={activeView === "registration" ? "workspace-nav-item workspace-nav-item-active" : "workspace-nav-item"}
            onClick={() => {
              setRequestedRegno(null);
              setActiveView("registration");
            }}
            type="button"
          >Register</button>
          <button
            aria-pressed={activeView === "metadata"}
            className={activeView === "metadata" ? "workspace-nav-item workspace-nav-item-active" : "workspace-nav-item"}
            onClick={() => setActiveView("metadata")}
            type="button"
          >Metadata</button>
          <button
            aria-pressed={activeView === "bulk"}
            className={activeView === "bulk" ? "workspace-nav-item workspace-nav-item-active" : "workspace-nav-item"}
            onClick={() => setActiveView("bulk")}
            type="button"
          >Bulk import</button>
        </nav>
        <div className="topbar-actions">
          <span className="database-badge"><Database size={14} />{session.database}</span>
          <span className="user-badge"><UserRound size={15} />{session.username}</span>
          <button className="icon-button logout-button" onClick={() => onLogout()} title="Sign out" type="button">
            <LogOut size={17} />
            <span>Sign out</span>
          </button>
        </div>
      </header>

      <section className="workspace-heading">
        <div>
          <span className="eyebrow">CHEMICAL REGISTRY / {session.database.toUpperCase()}</span>
          <h1>{activeView === "search" ? "Compound search" : activeView === "metadata" ? "Metadata" : activeView === "bulk" ? "Bulk import" : requestedRegno ? "Edit registration" : "Register compound"}</h1>
        </div>
        <div className="heading-mark" aria-hidden="true"><span>CR</span><i /></div>
      </section>

      {activeView === "search" ? <div className="workspace-grid">
        <aside className="search-panel">
          <div className="search-panel-heading">
            <div className="section-icon"><ScanSearch size={18} /></div>
            <div>
              <span className="eyebrow">RECORD INDEX</span>
              <h2>Find a compound</h2>
            </div>
          </div>

          <div className="search-modes" role="group" aria-label="Search field">
            {searchModes.map((item) => (
              <button
                aria-pressed={mode === item.value}
                className={mode === item.value ? "mode-button mode-button-active" : "mode-button"}
                key={item.value}
                onClick={() => {
                  setMode(item.value);
                  setSearchError("");
                }}
                type="button"
              >
                {item.label}
              </button>
            ))}
          </div>

          <form className="search-form" onSubmit={runSearch}>
            <label className="form-field" htmlFor="record-search">
              <span>{searchModes.find((item) => item.value === mode)?.label}</span>
              <span className="input-wrap search-input-wrap">
                <Search size={17} aria-hidden="true" />
                <input
                  autoComplete={mode === "submitter" ? "off" : "on"}
                  id="record-search"
                  list={mode === "submitter" ? "chemist-options" : undefined}
                  onChange={(event) => setSearchValue(event.target.value)}
                  placeholder={searchPlaceholders[mode]}
                  type="search"
                  value={searchValue}
                />
                {mode === "submitter" && (
                  <datalist id="chemist-options">
                    {chemists.map((chemist) => <option key={chemist} value={chemist} />)}
                  </datalist>
                )}
              </span>
            </label>
            {chemistError && mode === "submitter" && <p className="inline-note">{chemistError}</p>}
            {searchError && <p className="inline-error" role="alert">{searchError}</p>}
            <button className="primary-button search-submit" disabled={searchBusy} type="submit">
              {searchBusy ? <LoaderCircle className="spin" size={17} /> : <Search size={17} />}
              <span>{searchBusy ? "Searching" : "Search records"}</span>
            </button>
          </form>

          <div className="results-heading">
            <span>Results</span>
            <span className="result-count">{results.length.toString().padStart(2, "0")}</span>
          </div>

          {!searched && (
            <div className="results-empty">
              <span className="empty-mark"><Search size={18} /></span>
              <p>Search results will appear here.</p>
            </div>
          )}
          {searched && results.length === 0 && !searchBusy && !searchError && (
            <div className="results-empty">
              <span className="empty-mark"><ScanSearch size={18} /></span>
              <p>No matching registrations.</p>
            </div>
          )}
          {searchBusy && <div className="results-loading"><LoaderCircle className="spin" size={18} />Searching registry</div>}
          {results.length > 0 && (
            <ol className="result-list">
              {results.map((regno, index) => (
                <li key={`${regno}-${index}`}>
                  <button
                    aria-current={selectedRegno === regno ? "true" : undefined}
                    className={selectedRegno === regno ? "result-item result-item-active" : "result-item"}
                    onClick={() => setSelectedRegno(regno)}
                    type="button"
                  >
                    <span className="result-sequence">{String(index + 1).padStart(2, "0")}</span>
                    <span className="result-regno">{regno}</span>
                    <ChevronRight size={16} />
                  </button>
                </li>
              ))}
            </ol>
          )}
          <div className="search-panel-foot">
            <span>CHEMREG / INDEX</span>
            <Activity size={14} />
          </div>
        </aside>

        <section className="record-panel" aria-live="polite">
          {!selectedRegno && (
            <div className="record-placeholder">
              <div className="placeholder-symbol"><Atom size={28} strokeWidth={1.35} /></div>
              <span className="eyebrow">REGISTRATION RECORD</span>
              <h2>Record details</h2>
              <p>Select a registration from the results list.</p>
            </div>
          )}

          {selectedRegno && (
            <>
              <div className="record-toolbar">
                <div className="record-position">
                  <span className="eyebrow">REGISTRATION</span>
                  <span className="record-position-count">
                    {String(selectedIndex + 1).padStart(2, "0")} <i>/</i> {String(results.length).padStart(2, "0")}
                  </span>
                </div>
                <div className="record-nav">
                  <button
                    aria-label="Previous registration"
                    className="icon-button"
                    disabled={selectedIndex <= 0}
                    onClick={() => moveRecord(-1)}
                    title="Previous registration"
                    type="button"
                  ><ChevronLeft size={18} /></button>
                  <button
                    aria-label="Next registration"
                    className="icon-button"
                    disabled={selectedIndex < 0 || selectedIndex >= results.length - 1}
                    onClick={() => moveRecord(1)}
                    title="Next registration"
                    type="button"
                  ><ChevronRight size={18} /></button>
                </div>
              </div>

              {recordError && (
                <div className="record-error" role="alert">
                  <CircleAlert size={18} />
                  <span>{recordError}</span>
                  <button className="icon-button" onClick={() => setRecordRevision((value) => value + 1)} title="Retry loading record" type="button">
                    <RotateCcw size={16} />
                  </button>
                </div>
              )}

              {recordBusy && (
                <div className="record-loading"><LoaderCircle className="spin" size={20} /><span>Loading registration</span></div>
              )}

              {record && !recordBusy && (
                <RecordDetails
                  record={record}
                  imageUrl={molImageUrl}
                  imageFailed={molImageFailed}
                  onImageError={() => setMolImageFailed(true)}
                  onEdit={() => {
                    setRequestedRegno(record.regno);
                    setActiveView("registration");
                  }}
                />
              )}
            </>
          )}
        </section>
      </div> : activeView === "registration" ? (
        <RegistrationScreen
          key={requestedRegno ?? "new-registration"}
          requestedRegno={requestedRegno}
          session={session}
          onCancel={() => {
            setRequestedRegno(null);
            setActiveView("search");
          }}
          onComplete={(regno) => {
            setRequestedRegno(null);
            setResults([regno]);
            setSelectedRegno(regno);
            setSearched(true);
            setActiveView("search");
          }}
          onSessionExpired={expireSession}
        />
      ) : activeView === "metadata" ? (
        <MetadataScreen session={session} onSessionExpired={expireSession} />
      ) : (
        <BulkImportScreen session={session} onSessionExpired={expireSession} />
      )}
    </main>
  );
}

interface RecordDetailsProps {
  record: CompoundRecord;
  imageUrl: string;
  imageFailed: boolean;
  onImageError: () => void;
  onEdit: () => void;
}

function RecordDetails({ record, imageUrl, imageFailed, onImageError, onEdit }: RecordDetailsProps) {
  return (
    <div className="record-content">
      <div className="record-title-row">
        <div>
          <span className="eyebrow">CHEMICAL REGISTRATION</span>
          <h2>{record.compoundId || `Registration ${record.regno}`}</h2>
          <span className="record-id">REG NO. {record.regno}</span>
        </div>
        <div className="record-title-actions">
          <button className="secondary-button edit-record-button" onClick={onEdit} type="button">
            <Pencil size={15} /><span>Edit</span>
          </button>
          <span className="record-state">RECORD VIEW</span>
        </div>
      </div>

      <div className="record-overview-grid">
        <section className="structure-panel" aria-label="Molecule structure">
          <div className="section-kicker">
            <span>STRUCTURE</span>
            <span className="structure-status"><i />2D VIEW</span>
          </div>
          <div className="structure-canvas">
            {imageUrl && !imageFailed ? (
              <img alt={`Molecular structure for registration ${record.regno}`} onError={onImageError} src={imageUrl} />
            ) : (
              <div className="structure-empty">
                <Atom size={31} strokeWidth={1.25} />
                <span>{imageFailed ? "Structure image unavailable" : "Structure preview"}</span>
              </div>
            )}
          </div>
        </section>

        <section className="properties-panel" aria-label="Molecular properties">
          <div className="section-kicker"><span>MOLECULAR PROPERTIES</span><span className="property-mark">01</span></div>
          <div className="property-block property-formula">
            <span>MOLECULAR FORMULA</span>
            <strong>{record.formula || "-"}</strong>
          </div>
          <div className="property-pair">
            <div className="property-block">
              <span>AVERAGE MASS</span>
              <strong>{record.averageMass || "-"}</strong>
            </div>
            <div className="property-block">
              <span>MONOISOTOPIC MASS</span>
              <strong>{record.monoisotopicMass || "-"}</strong>
            </div>
          </div>
          <div className="property-block property-chns">
            <span>CHNS CALCULATION</span>
            <strong>{record.chns || "-"}</strong>
          </div>
        </section>
      </div>

      <section className="metadata-section">
        <div className="metadata-heading">
          <div>
            <span className="eyebrow">RECORD DATA</span>
            <h3>Registration metadata</h3>
          </div>
          <span className="metadata-date">{record.date || "DATE NOT AVAILABLE"}</span>
        </div>
        <dl className="metadata-grid">
          <MetadataField label="External ID" value={record.externalId} />
          <MetadataField label="External batch" value={record.externalBatch} />
          <MetadataField label="ELN batch" value={record.jpage} mono />
          <MetadataField label="Submitter" value={record.chemist} />
          <MetadataField label="Project" value={record.project} />
          <MetadataField label="Compound type" value={record.compoundType} />
          <MetadataField label="Product type" value={record.product} />
          <MetadataField label="Library" value={record.library} />
          <MetadataField label="Analytical solvent" value={record.solvent} />
          <MetadataField label="IP rights" value={record.ipRights} />
        </dl>
      </section>
    </div>
  );
}

function MetadataField({ label, value, mono = false }: { label: string; value: string; mono?: boolean }) {
  return (
    <div className="metadata-field">
      <dt>{label}</dt>
      <dd className={mono ? "mono-value" : undefined}>{value || "-"}</dd>
    </div>
  );
}