import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

const proxyTarget = process.env.CHEMREG_API_TARGET ?? "https://esox3.scilifelab.se/chemreg/";
const proxy = { target: proxyTarget, changeOrigin: true };

export default defineConfig({
  plugins: [react()],
  server: {
    host: "0.0.0.0",
    port: 8082,
    strictPort: true,
    proxy: {
      "/login": proxy,
      "/getDatabase": proxy,
      "/api": proxy,
      "/mols": proxy
    }
  },
  preview: {
    host: "0.0.0.0",
    port: 8082,
    strictPort: true
  },
  build: {
    outDir: "build",
    emptyOutDir: true
  }
});