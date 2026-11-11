import { afterEach, describe, expect, it, vi } from "vitest";
import {
  login,
  createDraftRegno,
  createLibrary,
  parseStringList,
  parseStringValue,
  postChemRegMol,
  searchCompounds,
  updateCompoundValue,
  uploadMolfile
} from "./api";

afterEach(() => {
  vi.unstubAllGlobals();
});

describe("legacy API response parsing", () => {
  it("unwraps database and search rows", () => {
    expect(parseStringList([["Live"], ["Test"]])).toEqual(["Live", "Test"]);
    expect(parseStringList([[123], [null], ["None"]])).toEqual(["123"]);
  });

  it("unwraps a single value returned as a nested row", () => {
    expect(parseStringValue([["CHEM-42"]])).toBe("CHEM-42");
    expect(parseStringValue([])).toBe("");
  });

  it("keeps the exact login response text for the legacy token header", async () => {
    const responseText = '{"token": "signed-token", "database": "Test"}';
    const fetchMock = vi.fn().mockResolvedValue(
      new Response(responseText, { status: 200 })
    );
    vi.stubGlobal("fetch", fetchMock);

    const session = await login("chemist", "secret", "Test");

    expect(session.tokenHeader).toBe(responseText);
    expect(fetchMock).toHaveBeenCalledWith(
      "/login",
      expect.objectContaining({ method: "POST" })
    );
  });

  it("sends the raw login response as the protected API token header", async () => {
    const tokenHeader = '{"token": "signed-token", "database": "Test"}';
    const fetchMock = vi.fn().mockResolvedValue(
      new Response('[["70001"]]', { status: 200 })
    );
    vi.stubGlobal("fetch", fetchMock);

    await searchCompounds("regno", "70001", {
      username: "chemist",
      database: "Test",
      tokenHeader
    });

    const [, options] = fetchMock.mock.calls[0] as [string, RequestInit];
    expect(new Headers(options.headers).get("token")).toBe(tokenHeader);
  });

  it("searches compound IDs with the same column key as the desktop client", async () => {
    const fetchMock = vi.fn().mockResolvedValue(
      new Response('[["70001"]]', { status: 200 })
    );
    vi.stubGlobal("fetch", fetchMock);

    const hits = await searchCompounds("compoundId", "CBK621542", {
      username: "chemist",
      database: "Test",
      tokenHeader: '{"token": "signed-token", "database": "Test"}'
    });

    expect(fetchMock.mock.calls[0]?.[0]).toBe(
      "/api/search?column=COMPOUND_ID&value=CBK621542"
    );
    expect(hits).toEqual(["70001"]);
  });

  it("creates a draft through next-regno then create-regno", async () => {
    const fetchMock = vi.fn()
      .mockResolvedValueOnce(new Response('"70001"', { status: 200 }))
      .mockResolvedValueOnce(new Response("", { status: 200 }));
    vi.stubGlobal("fetch", fetchMock);

    const regno = await createDraftRegno({
      username: "chemist",
      database: "Test",
      tokenHeader: '{"token": "signed-token", "database": "Test"}'
    });

    expect(regno).toBe("70001");
    expect(fetchMock.mock.calls[0]?.[0]).toBe("/api/getNextRegno");
    expect(fetchMock.mock.calls[1]?.[0]).toBe("/api/createRegno?regno=70001");
    expect(fetchMock.mock.calls[1]?.[1]).toEqual(expect.objectContaining({ method: "PUT" }));
  });

  it("updates fields with the legacy PUT query contract", async () => {
    const fetchMock = vi.fn().mockResolvedValue(new Response("", { status: 200 }));
    vi.stubGlobal("fetch", fetchMock);

    await updateCompoundValue("70001", "CHEMIST", "Ada Chemist", {
      username: "chemist",
      database: "Test",
      tokenHeader: "legacy-token"
    });

    expect(fetchMock.mock.calls[0]?.[0]).toBe(
      "/api/update?column=CHEMIST&value=Ada+Chemist&regno=70001"
    );
    expect(fetchMock.mock.calls[0]?.[1]).toEqual(expect.objectContaining({ method: "PUT" }));
  });

  it("sends the regno as a multipart file part for molfile upload", async () => {
    const fetchMock = vi.fn().mockResolvedValue(new Response("", { status: 200 }));
    vi.stubGlobal("fetch", fetchMock);
    const file = new File(["molfile content"], "sample.mol", { type: "chemical/x-mdl-molfile" });

    await uploadMolfile("70001", file, {
      username: "chemist",
      database: "Test",
      tokenHeader: "legacy-token"
    });

    const options = fetchMock.mock.calls[0]?.[1] as RequestInit;
    const body = options.body as FormData;
    expect(fetchMock.mock.calls[0]?.[0]).toBe("/api/loadMolfile");
    const uploadedFile = body.get("file") as File;
    expect(uploadedFile.name).toBe("sample.mol");
    expect(await uploadedFile.text()).toBe("molfile content");
    expect(await (body.get("regno") as File).text()).toBe("70001");
  });

  it("uses the server's capitalized restricted-library flag", async () => {
    const fetchMock = vi.fn().mockResolvedValue(new Response("", { status: 200 }));
    vi.stubGlobal("fetch", fetchMock);

    await createLibrary("Library A", "Supplier A", false, {
      username: "chemist",
      database: "Test",
      tokenHeader: "legacy-token"
    });

    expect(fetchMock.mock.calls[0]?.[0]).toContain("restricted=False");
  });

  it("posts molecule fields as form data and keeps row-level HTTP failures inspectable", async () => {
    const fetchMock = vi.fn().mockResolvedValue(
      new Response("Bad molecule", { status: 500 })
    );
    vi.stubGlobal("fetch", fetchMock);

    const result = await postChemRegMol({ smiles: "CCO", jpage: "ABC123001" }, {
      username: "chemist",
      database: "Test",
      tokenHeader: "legacy-token"
    });

    const options = fetchMock.mock.calls[0]?.[1] as RequestInit;
    expect(fetchMock.mock.calls[0]?.[0]).toBe("/api/chemRegAddMol");
    expect((options.body as URLSearchParams).get("smiles")).toBe("CCO");
    expect(result).toEqual({ ok: false, status: 500, body: "Bad molecule" });
  });
});