export interface Session {
  username: string;
  database: string;
  tokenHeader: string;
}

export type SearchMode = "regno" | "submitter" | "compoundId" | "batch";

export interface CompoundRecord {
  regno: string;
  date: string;
  compoundId: string;
  externalId: string;
  externalBatch: string;
  jpage: string;
  chemist: string;
  project: string;
  compoundType: string;
  product: string;
  library: string;
  solvent: string;
  ipRights: string;
  formula: string;
  averageMass: string;
  monoisotopicMass: string;
  chns: string;
}

export interface RegistrationRecord extends CompoundRecord {
  purity: string;
  chromText: string;
  chromPurity: string;
  nmrText: string;
  nmrPurity: string;
  msText: string;
  msPurity: string;
  comments: string;
}

export interface ApiResult {
  ok: boolean;
  status: number;
  body: string;
}

export class ApiError extends Error {
  constructor(
    message: string,
    public readonly status: number
  ) {
    super(message);
    this.name = "ApiError";
  }
}

const searchColumns: Record<SearchMode, string> = {
  regno: "regno",
  submitter: "CHEMIST",
  compoundId: "COMPOUND_ID",
  batch: "JPAGE"
};

const recordColumns = [
  ["date", "rdate"],
  ["compoundId", "compound_id"],
  ["externalId", "EXTERNAL_ID"],
  ["externalBatch", "SUPPLIER_BATCH"],
  ["jpage", "JPAGE"],
  ["chemist", "chemist"],
  ["project", "project"],
  ["compoundType", "compound_type"],
  ["product", "product"],
  ["library", "library_id"],
  ["solvent", "solvent"],
  ["ipRights", "ip_rights"],
  ["formula", "C_MF"],
  ["averageMass", "C_MW"],
  ["monoisotopicMass", "C_MONOISO"],
  ["chns", "C_CHNS"]
] as const satisfies ReadonlyArray<readonly [keyof CompoundRecord, string]>;

const registrationColumns = [
  ["purity", "PURITY"],
  ["chromText", "CHROM_TEXT"],
  ["chromPurity", "CHROM_PURITY"],
  ["nmrText", "NMR_TEXT"],
  ["nmrPurity", "NMR_PURITY"],
  ["msText", "MS_TEXT"],
  ["msPurity", "MS_PURITY"],
  ["comments", "COMMENTS"]
] as const satisfies ReadonlyArray<readonly [keyof RegistrationRecord, string]>;

function parseBody(text: string): unknown {
  if (!text.trim()) return null;
  try {
    return JSON.parse(text) as unknown;
  } catch {
    return text;
  }
}

function responseError(text: string, status: number): ApiError {
  const body = parseBody(text);
  if (typeof body === "string" && body.trim()) {
    return new ApiError(body, status);
  }
  if (typeof body === "object" && body !== null && "message" in body) {
    const message = (body as { message?: unknown }).message;
    if (typeof message === "string" && message.trim()) {
      return new ApiError(message, status);
    }
  }
  return new ApiError(`Request failed (${status}).`, status);
}

async function request(
  path: string,
  session?: Session,
  init: RequestInit = {}
): Promise<unknown> {
  const headers = new Headers(init.headers);
  if (session) headers.set("token", session.tokenHeader);
  const response = await fetch(path, { ...init, headers });
  const text = await response.text();
  if (!response.ok) throw responseError(text, response.status);
  return parseBody(text);
}

async function postForm(
  path: string,
  values: Record<string, string>,
  session: Session
): Promise<ApiResult> {
  const headers = new Headers({ "Content-Type": "application/x-www-form-urlencoded" });
  headers.set("token", session.tokenHeader);
  const response = await fetch(path, {
    method: "POST",
    headers,
    body: new URLSearchParams(values)
  });
  return { ok: response.ok, status: response.status, body: await response.text() };
}

export function parseStringList(value: unknown): string[] {
  let rows = value;
  if (typeof rows === "string") {
    const text = rows;
    try {
      rows = JSON.parse(text) as unknown;
    } catch {
      return text.trim() ? [text] : [];
    }
  }
  if (!Array.isArray(rows)) return [];

  return rows.flatMap((row) => {
    const cell = Array.isArray(row) ? row[0] : row;
    if (cell === null || cell === undefined) return [];
    const text = String(cell).trim();
    return text && text.toLowerCase() !== "none" ? [text] : [];
  });
}

export function parseStringValue(value: unknown): string {
  if (Array.isArray(value)) {
    return value.length ? parseStringValue(value[0]) : "";
  }
  if (value === null || value === undefined) return "";
  if (typeof value === "object") {
    const firstValue = Object.values(value)[0];
    return firstValue === undefined ? "" : parseStringValue(firstValue);
  }
  return String(value).trim();
}

export async function getDatabases(): Promise<string[]> {
  return parseStringList(await request("/getDatabase"));
}

export async function login(
  username: string,
  password: string,
  database: string
): Promise<Session> {
  const response = await fetch("/login", {
    method: "POST",
    headers: { "Content-Type": "application/x-www-form-urlencoded" },
    body: new URLSearchParams({ username, password, database })
  });
  const tokenHeader = await response.text();
  if (!response.ok) throw responseError(tokenHeader, response.status);

  const payload = parseBody(tokenHeader);
  if (
    typeof payload !== "object" ||
    payload === null ||
    !("token" in payload) ||
    typeof (payload as { token?: unknown }).token !== "string"
  ) {
    throw new ApiError("The server returned an invalid login response.", response.status);
  }
  return { username, database, tokenHeader };
}

export async function getColumnValues(
  column: string,
  session: Session
): Promise<string[]> {
  const params = new URLSearchParams({ column });
  return parseStringList(
    await request(`/api/getColComboData?${params.toString()}`, session)
  );
}

export async function searchCompounds(
  mode: SearchMode,
  value: string,
  session: Session
): Promise<string[]> {
  const params = new URLSearchParams({ column: searchColumns[mode], value });
  return parseStringList(await request(`/api/search?${params.toString()}`, session));
}

export async function getCompoundRecord(
  regno: string,
  session: Session
): Promise<CompoundRecord> {
  return { regno, ...await getRecordFields(regno, session, recordColumns) } as CompoundRecord;
}

export async function getRegistrationRecord(
  regno: string,
  session: Session
): Promise<RegistrationRecord> {
  return {
    regno,
    ...await getRecordFields(regno, session, [...recordColumns, ...registrationColumns])
  } as RegistrationRecord;
}

async function getRecordFields(
  regno: string,
  session: Session,
  columns: ReadonlyArray<readonly [string, string]>
): Promise<Record<string, string>> {
  const fields = await Promise.all(
    columns.map(async ([key, column]) => {
      const params = new URLSearchParams({ column, regno });
      const value = await request(
        `/api/getTextColumn?${params.toString()}`,
        session
      );
      return [key, parseStringValue(value)] as const;
    })
  );
  return Object.fromEntries(fields);
}

export async function createMolImageUrl(
  regno: string,
  session: Session
): Promise<string> {
  const params = new URLSearchParams({ regno });
  await request(`/api/createMolImage?${params.toString()}`, session);
  return `/mols/${encodeURIComponent(regno)}.png?v=${Date.now()}`;
}

export async function createDraftRegno(session: Session): Promise<string> {
  const regno = parseStringValue(await request("/api/getNextRegno", session));
  if (!regno) throw new ApiError("The server did not return a registration number.", 500);
  const params = new URLSearchParams({ regno });
  await request(`/api/createRegno?${params.toString()}`, session, { method: "PUT" });
  return regno;
}

export async function deleteDraftRegno(regno: string, session: Session): Promise<void> {
  const params = new URLSearchParams({ regno });
  await request(`/api/deleteRegno?${params.toString()}`, session, { method: "PUT" });
}

export async function updateCompoundValue(
  regno: string,
  column: string,
  value: string,
  session: Session
): Promise<void> {
  const params = new URLSearchParams({ column, value, regno });
  await request(`/api/update?${params.toString()}`, session, { method: "PUT" });
}

export async function updateRegistrationBatch(
  regno: string,
  batch: string,
  session: Session
): Promise<void> {
  const params = new URLSearchParams({ regno, batch });
  await request(`/api/updateRegnoBatch?${params.toString()}`, session, { method: "PUT" });
}

export async function uploadMolfile(
  regno: string,
  file: File,
  session: Session
): Promise<void> {
  const formData = new FormData();
  formData.append("file", file, file.name);
  formData.append("regno", new Blob([regno], { type: "text/plain" }), "regno.txt");
  await request("/api/loadMolfile", session, { method: "POST", body: formData });
}

export async function getMolfile(regno: string, session: Session): Promise<string> {
  const params = new URLSearchParams({ regno });
  return parseStringValue(await request(`/api/getMolfile?${params.toString()}`, session));
}

export async function registerCompound(regno: string, session: Session): Promise<void> {
  const params = new URLSearchParams({ regno });
  await request(`/api/bcpvsRegCompound?${params.toString()}`, session, { method: "PUT" });
}

export async function getCanonicalSmiles(
  smiles: string,
  session: Session
): Promise<{ exists: boolean; smiles: string }> {
  const params = new URLSearchParams({ smiles });
  const response = await request(`/api/getCanonicSmiles?${params.toString()}`, session);
  const values = Array.isArray(response) ? response : [];
  const exists = parseStringValue(values[0]) !== "False";
  return { exists, smiles: parseStringValue(values[1]) };
}

export async function createSalt(smiles: string, session: Session): Promise<void> {
  const params = new URLSearchParams({ smiles });
  await request(`/api/createSalt?${params.toString()}`, session, { method: "PUT" });
}

export async function createSupplier(supplier: string, session: Session): Promise<void> {
  const params = new URLSearchParams({ supplier });
  await request(`/api/createSupplier?${params.toString()}`, session, { method: "PUT" });
}

export async function createLibrary(
  libraryName: string,
  supplier: string,
  restricted: boolean,
  session: Session
): Promise<void> {
  const params = new URLSearchParams({
    library_name: libraryName,
    supplier,
    restricted: restricted ? "True" : "False"
  });
  await request(`/api/createLibrary?${params.toString()}`, session, { method: "PUT" });
}

export async function getLastBatchOfEln(eln: string, session: Session): Promise<number> {
  const params = new URLSearchParams({ eln });
  const value = Number(parseStringValue(
    await request(`/api/getLastBatchFromEln?${params.toString()}`, session)
  ));
  if (!Number.isInteger(value) || value < 0) {
    throw new ApiError(`The server returned an invalid batch count for ${eln}.`, 500);
  }
  return value;
}

export async function getLastCtrlCompound(session: Session): Promise<number> {
  const value = Number(parseStringValue(await request("/api/getLastCtrlComp?eln=CTRL", session)));
  if (!Number.isInteger(value) || value < 0) {
    throw new ApiError("The server returned an invalid control-compound count.", 500);
  }
  return value;
}

export async function getNextSdfSequence(session: Session): Promise<string> {
  const value = parseStringValue(await request("/api/getNextSdfSequence", session));
  if (!value) throw new ApiError("The server did not return an SDF sequence number.", 500);
  return value;
}

export async function getRegnosFromSdfSequence(
  sequence: string,
  session: Session
): Promise<string[]> {
  const params = new URLSearchParams({ sdfile_sequence: sequence });
  return parseStringList(
    await request(`/api/getRegnosFromSequence?${params.toString()}`, session)
  );
}

export async function postChemRegMol(
  values: Record<string, string>,
  session: Session
): Promise<ApiResult> {
  return postForm("/api/chemRegAddMol", values, session);
}

export async function postNostructMol(
  values: Record<string, string>,
  session: Session
): Promise<ApiResult> {
  return postForm("/api/addNostructMol", values, session);
}

export async function postCtrlMol(
  values: Record<string, string>,
  session: Session
): Promise<ApiResult> {
  return postForm("/api/addCtrlMol", values, session);
}