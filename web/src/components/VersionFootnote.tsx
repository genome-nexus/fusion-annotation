/** Small build-provenance footnote — shows the released package version and
 * the exact commit the running bundle was built from. Baked in at build
 * time via VITE_APP_VERSION / VITE_COMMIT_SHA (see deploy-on-release.yml /
 * deploy-api-web.yml); both are undefined for a plain local `npm run dev`,
 * where we fall back to "dev" and hide the commit link. Handy for quickly
 * confirming what's actually live without digging through Actions runs —
 * see the GitHub Pages stale-build incident this was added after. */
export function VersionFootnote() {
  const version = import.meta.env.VITE_APP_VERSION;
  const commit = import.meta.env.VITE_COMMIT_SHA;

  return (
    <footer className="version-footnote">
      {version ? `v${version}` : "dev build"}
      {commit && (
        <>
          {" · "}
          <a
            href={`https://github.com/genome-nexus/fusion-annotation/commit/${commit}`}
            target="_blank"
            rel="noreferrer"
          >
            {commit.slice(0, 7)}
          </a>
        </>
      )}
    </footer>
  );
}
