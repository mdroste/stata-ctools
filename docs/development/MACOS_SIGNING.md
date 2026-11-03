# macOS plugin signing and notarization

`Build Plugins` signs and notarizes both macOS plugins before packaging. Builds
on `main`, `dev`, and `develop`, and releases with `v`-prefixed tags, use the
`macos-signing` GitHub environment. Pull requests build candidates without Apple
credentials or Developer ID signatures (ARM binaries may have linker-generated
ad hoc signatures). Other manually selected refs cannot produce a distribution.

## Apple credentials

An Apple Developer account or a certificate installed on your Mac does not make
its signing identity available to GitHub's runners. The workflow needs an
exportable **Developer ID Application** certificate together with its private
key, plus credentials for Apple's notary service. Use an active Apple Developer
Program membership and the team that owns that certificate.

1. Open Keychain Access on the Mac that holds the certificate's private key.
   Under **My Certificates**, locate **Developer ID Application: Your Name
   (TEAMID)**. Expand it and confirm that a private key appears. Export the
   identity as `DeveloperIDApplication.p12`, including the private key, and set a
   strong export password. A downloaded `.cer` alone cannot sign code. If the
   private key is missing, obtain the original identity or create a new
   Developer ID Application certificate through your team's Account Holder.
   Do not use Developer ID Installer, Apple Development, or Mac App Distribution.
2. Encode the exported identity for the GitHub secret:
   ```bash
   base64 -i /absolute/path/to/DeveloperIDApplication.p12 | pbcopy
   ```
3. Find the ten-character **Team ID** in your Apple Developer account's
   membership details. It must match the certificate's team.
4. At [account.apple.com](https://account.apple.com/), choose **Sign-In and
   Security → App-Specific Passwords**, and generate a password for ctools CI.
   Use the Apple Account associated with that developer team. Two-factor
   authentication is required. This is separate from both your ordinary Apple
   Account password and the `.p12` export password. Changing your ordinary
   account password revokes existing app-specific passwords.

See Apple's [Developer ID signing instructions](https://developer.apple.com/documentation/xcode/creating-distribution-signed-code-for-the-mac),
[notarization workflow](https://developer.apple.com/documentation/security/customizing-the-notarization-workflow),
and [app-specific password instructions](https://support.apple.com/en-us/102654),
and GitHub's [certificate import guidance](https://docs.github.com/en/actions/how-tos/deploy/deploy-to-third-party-platforms/sign-xcode-applications).

## GitHub environment

In repository **Settings → Environments**, create `macos-signing`. Add:

| Kind | Name | Value |
|---|---|---|
| Secret | `APPLE_CERTIFICATE_P12_BASE64` | Base64 of the exported Developer ID Application certificate and private key |
| Secret | `APPLE_CERTIFICATE_PASSWORD` | Password used to encrypt the `.p12` export |
| Secret | `APPLE_ID` | Apple Account email used for notarization |
| Secret | `APPLE_APP_SPECIFIC_PASSWORD` | App-specific password generated for this workflow |
| Variable | `APPLE_TEAM_ID` | Ten-character Developer Team ID matching the certificate |

Restrict environment deployment branches/tags to `main`, `dev`, `develop`, and
release tags matching `v*`. Only trusted maintainers should be able to change
code on those refs. Pull requests never enter the signing job. No provisioning
profile, Developer ID Installer certificate, App Store listing, or saved
keychain password is needed for these plugins; CI generates an ephemeral
keychain password. Store credentials directly in GitHub, never in the repository.

## Build and publication

The compiler jobs upload `macos-arm-unsigned` and `macos-intel-unsigned`. A separate
job downloads those artifacts from the same workflow run and invokes
`scripts/sign_macos.sh`. It imports the identity into a temporary keychain,
requires exactly one valid Developer ID Application identity for the specified
team, signs each plugin with a secure timestamp, and checks its signature and
team. These Mach-O bundles run inside Stata; they need no app entitlements or
Hardened Runtime flag of their own.

CI submits a ZIP containing both signed plugins with `notarytool`, waits up to
30 minutes, checks for `Accepted` in both the response and service log, and
verifies Apple's online notarization ticket for each binary. The response and
log are saved in the `macos-notarization` Actions artifact, including on failure
when available. A timeout or failed signature/ticket check blocks distribution.
A timeout does not cancel Apple's processing; its submission ID can be used to
inspect the request. The temporary keychain and certificate files are removed
on exit. Signing credentials are limited to this step on a disposable hosted
runner.

Only after verification does CI upload `ctools_mac_arm` and `ctools_mac_x86`.
Packaging downloads only the `ctools_*` artifacts, requires both macOS and
Windows signing jobs to succeed outside pull requests, and calculates checksums
from the final signed bytes. Missing credentials fail the job; there is no
unsigned fallback for branch or release distributions.

Once the workflow change is on `main` and credentials are configured, run
**Actions → Build Plugins → Run workflow** on `main` with **publish_plugins**
selected. All platform checks, including Windows signing, must pass. This updates
the binaries served by `net install`; pushing workflow changes alone does not
replace them. Creating a `v`-prefixed GitHub release builds and attaches a new
release archive; existing release downloads are not retroactively replaced.

## Verification and limitations

Apple issues tickets for standalone Mach-O files, but cannot staple tickets to
those files or directly to ZIP archives. These `.plugin` files are standalone
binaries, not directory bundles. Keep the signed bytes unchanged when copying
or packaging them. Gatekeeper may need internet access on first use to retrieve
the ticket. Fully offline first-use support would require a different delivery
container, such as a stapled installer package or disk image.

For an installed plugin, inspect and verify it in Terminal:

```bash
codesign -dvvv /absolute/path/to/ctools_mac_x86.plugin
codesign --verify --strict --verbose=2 /absolute/path/to/ctools_mac_x86.plugin
codesign --verify --verbose=2 --check-notarization -R='notarized' /absolute/path/to/ctools_mac_x86.plugin
```

Use `ctools_mac_arm.plugin` for Apple Silicon. Expect a Developer ID Application
authority, the configured Team Identifier, and a timestamp. Test the downloaded
package in Stata on clean Intel and Apple Silicon Macs with Gatekeeper enabled
and internet access. These workflow checks do not exercise licensed Stata or
prove runtime compatibility on every supported OS. A normal signature check
alone does not establish notarization, and `spctl --assess --type execute` is
intended for applications rather than these plugin files.

For import failures, confirm the `.p12` contains the private key and the export
password is correct. For identity failures, check certificate type, expiration,
and Team ID. For authentication failures, check team membership and the
app-specific password. For rejection, inspect `notarization-log.json`; for ticket
lookup failures after acceptance, check Apple service availability and rerun.
Renew the certificate and replace the two certificate secrets before it expires.
