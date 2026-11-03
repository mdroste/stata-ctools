#!/bin/bash
# Run only on an ephemeral macOS runner; never enable shell tracing here.
set -euo pipefail
umask 077

for name in APPLE_CERTIFICATE_P12_BASE64 APPLE_CERTIFICATE_PASSWORD APPLE_ID \
            APPLE_APP_SPECIFIC_PASSWORD APPLE_TEAM_ID; do
    if [[ -z "${!name:-}" ]]; then
        echo "Set $name in the macos-signing GitHub environment; see docs/development/MACOS_SIGNING.md." >&2
        exit 1
    fi
done
if [[ ! "$APPLE_TEAM_ID" =~ ^[A-Z0-9]{10}$ ]]; then
    echo 'APPLE_TEAM_ID must be the 10-character Apple Developer Team ID.' >&2
    exit 1
fi

plugins="${1:?Usage: sign_macos.sh PLUGIN_DIRECTORY RESULT_DIRECTORY}"
results="${2:?Usage: sign_macos.sh PLUGIN_DIRECTORY RESULT_DIRECTORY}"
for arch in arm x86; do
    test -f "$plugins/ctools_mac_$arch.plugin"
done
mkdir -p "$results"
work=$(mktemp -d "${RUNNER_TEMP:?}/ctools-signing.XXXXXX")
keychain="$work/signing.keychain-db"
cleanup() {
    # Preserve the original failure, but never report success if cleanup fails.
    result=$?
    trap - EXIT
    if [[ -f "$keychain" ]]; then
        security delete-keychain "$keychain" || result=1
    fi
    rm -rf "$work" || result=1
    exit "$result"
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

keychain_password=$(openssl rand -hex 32)
printf '%s' "$APPLE_CERTIFICATE_P12_BASE64" | base64 --decode > "$work/certificate.p12"
security create-keychain -p "$keychain_password" "$keychain"
security set-keychain-settings -lut 3600 "$keychain"
security unlock-keychain -p "$keychain_password" "$keychain"
security import "$work/certificate.p12" -P "$APPLE_CERTIFICATE_PASSWORD" \
    -t cert -f pkcs12 -k "$keychain" -T /usr/bin/codesign -T /usr/bin/security
security set-key-partition-list -S apple-tool:,apple:,codesign: -s \
    -k "$keychain_password" "$keychain" > /dev/null
rm -f "$work/certificate.p12"

# Select exactly one usable Developer ID Application identity from this keychain.
security find-identity -v -p codesigning "$keychain" > "$work/identities.txt"
identity=$(python3 - "$work/identities.txt" <<'PY'
import os
import re
import sys
from pathlib import Path
identities = re.findall(
    r'\b([A-Fa-f0-9]{40}) "Developer ID Application: [^"\n]+ \('
    + re.escape(os.environ['APPLE_TEAM_ID']) + r'\)"',
    Path(sys.argv[1]).read_text(),
)
if len(identities) != 1:
    sys.exit('Export exactly one valid Developer ID Application certificate and its private key for APPLE_TEAM_ID.')
print(identities[0])
PY
)
requirement="anchor apple generic and certificate leaf[field.1.2.840.113635.100.6.1.13] exists and certificate leaf[subject.OU] = \"$APPLE_TEAM_ID\""
mkdir "$work/payload"
for arch in arm x86; do
    plugin="$plugins/ctools_mac_$arch.plugin"
    # MH_BUNDLE plugins inherit Stata's runtime policy; no app entitlements or
    # provisioning profile are needed. Hardened Runtime is for main executables.
    codesign --force --sign "$identity" --keychain "$keychain" --timestamp \
        --identifier "org.ctools.mac.$arch" "$plugin"
    codesign --verify --strict --verbose=2 "-R=$requirement" "$plugin"
    codesign --display --verbose=4 "$plugin" 2> "$work/signature.txt"
    if ! grep -q '^Timestamp=' "$work/signature.txt"; then
        echo "$plugin: missing secure timestamp" >&2
        exit 1
    fi
    cp "$plugin" "$work/payload/"
done

xcrun notarytool store-credentials ctools-notary --keychain "$keychain" \
    --apple-id "$APPLE_ID" --team-id "$APPLE_TEAM_ID" \
    --password "$APPLE_APP_SPECIFIC_PASSWORD"
ditto -c -k --keepParent "$work/payload" "$work/ctools-macos.zip"
# Keep the response and log even for rejection or timeout; neither permits upload.
submit_result=0
xcrun notarytool submit "$work/ctools-macos.zip" --keychain "$keychain" \
    --keychain-profile ctools-notary --wait --timeout 30m --output-format json \
    > "$results/submission.json" || submit_result=$?
submission_id=$(python3 - "$results/submission.json" <<'PY'
import json
import sys
import uuid
from pathlib import Path
response = json.loads(Path(sys.argv[1]).read_text())
print(uuid.UUID(response['id']))
PY
)
xcrun notarytool log "$submission_id" --keychain "$keychain" \
    --keychain-profile ctools-notary "$results/notarization-log.json"
python3 - "$results/submission.json" "$results/notarization-log.json" <<'PY'
import json
import sys
from pathlib import Path
for name in sys.argv[1:]:
    response = json.loads(Path(name).read_text())
    if response.get('status') != 'Accepted':
        sys.exit(f'Notarization was not accepted; inspect {name}.')
log = json.loads(Path(sys.argv[2]).read_text())
for issue in log.get('issues') or []:
    print(f"Notary service {issue.get('severity')}: {issue.get('message')}")
PY
if [[ "$submit_result" != 0 ]]; then
    exit "$submit_result"
fi
for arch in arm x86; do
    plugin="$plugins/ctools_mac_$arch.plugin"
    # Standalone Mach-O files cannot be stapled. Verify Apple's online ticket;
    # allow a short propagation interval after an Accepted submission.
    verified=false
    for attempt in 1 2 3 4 5 6; do
        if codesign --verify --strict --verbose=2 --check-notarization \
            "-R=($requirement) and notarized" "$plugin"; then
            verified=true
            break
        fi
        if [[ "$attempt" != 6 ]]; then sleep 10; fi
    done
    if [[ "$verified" != true ]]; then
        echo "$plugin: notarization ticket verification failed" >&2
        exit 1
    fi
    shasum -a 256 "$plugin"
done
if [[ -n "${GITHUB_STEP_SUMMARY:-}" ]]; then
    {
        echo 'macOS ARM and Intel plugins: Developer ID signatures, timestamps, and notarization verified.'
        echo "Apple Team ID: $APPLE_TEAM_ID"
        echo "Notarization submission: $submission_id"
        shasum -a 256 "$plugins/ctools_mac_arm.plugin" "$plugins/ctools_mac_x86.plugin"
    } >> "$GITHUB_STEP_SUMMARY"
fi
