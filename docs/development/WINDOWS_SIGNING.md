# Windows plugin signing

The `Build Plugins` workflow signs the Windows plugin with Azure Artifact Signing
before packaging. Branch builds on `main`, `dev`, and
`develop`, and releases with `v`-prefixed tags, use the `windows-signing` GitHub
environment. Pull requests produce unsigned candidate artifacts without Azure
access. Other manually selected refs do not produce a distributable package.

The Windows compiler uploads `windows-unsigned`. A separate Windows signing job
downloads that artifact from the same run, uses GitHub OIDC to authenticate, signs
and timestamps the PE DLL, and verifies Windows trust and timestamp presence.
Only then does it upload `ctools_windows`. Signing tools operate on a temporary
`.dll` filename; renaming it back to `.plugin` preserves the signed bytes. Package
checksums are calculated after signing. Missing configuration or failed signing
blocks packaging outside pull requests. Publication requires successful builds on all four platforms and the native and
platform dependency checks. GitHub Actions does not run licensed Stata tests.

## Azure setup

Create a Basic Artifact Signing account and complete public identity validation.
Create a **Public Trust** certificate profile, with **Program Type: None**.
The same profile can sign multiple projects under the verified publisher identity.

In **App registrations**, create a single-tenant application for this repository,
for example `ctools-github-signing`, with no redirect URI. Record its Application
(client) ID and Directory (tenant) ID. No client secret is needed.

Open that app's **Certificates & secrets > Federated credentials > Add credential**.
Select **GitHub actions deploying Azure resources** and enter:

| Field | Value |
|---|---|
| Organization | `mdroste` |
| Repository | `stata-ctools` |
| Entity type | `Environment` |
| Environment | `windows-signing` |
| Credential name | `ctools-windows-signing` |
| Issuer | `https://token.actions.githubusercontent.com` |
| Audience | `api://AzureADTokenExchange` |
| Subject identifier | `repo:mdroste@23185385/stata-ctools@1127448459:environment:windows-signing` |

The repository uses GitHub's immutable OIDC subject format. The numeric suffixes
are the owner and repository IDs. Keep them in Azure's generated subject. For a
new repository, enable immutable subjects in its OIDC settings and confirm the
`sub_claim_prefix` returned by `GET /repos/{owner}/{repo}/actions/oidc/customization/sub`
before configuring Azure. Append `:environment:windows-signing` to that prefix.

On the Artifact Signing account, open **Access control (IAM) > Add role assignment**.
Assign **Artifact Signing Certificate Profile Signer** to the application above
(search for its name under **User, group, or service principal**). The account scope
is sufficient; subscription-wide Contributor or Owner permissions are unnecessary.

See Microsoft's [account and profile setup](https://learn.microsoft.com/en-us/azure/artifact-signing/quickstart),
[GitHub federation instructions](https://learn.microsoft.com/en-us/entra/workload-id/workload-identity-federation-create-trust#github-actions),
and [signing roles](https://learn.microsoft.com/en-us/azure/artifact-signing/tutorial-assign-roles).

## GitHub setup

In repository **Settings > Environments**, create `windows-signing`. Under
**Deployment branches and tags**, select **Selected branches and tags**. Add
individual Branch rules for `main`, `dev`, and `develop`, and a Tag rule for `v*`.
Do not allow pull-request merge refs to use this environment.

Add these **environment secrets**:

| Secret | Value |
|---|---|
| `AZURE_CLIENT_ID` | Application (client) ID of `ctools-github-signing` |
| `AZURE_TENANT_ID` | Directory (tenant) ID of that application |
| `AZURE_SUBSCRIPTION_ID` | ID of the subscription containing the signing account |

Add these **environment variables**:

| Variable | Value |
|---|---|
| `AZURE_SIGNING_ACCOUNT` | Artifact Signing account name |
| `AZURE_SIGNING_PROFILE` | Public Trust certificate profile name |
| `AZURE_SIGNING_ENDPOINT` | Signing endpoint for the account's region; East US uses `https://eus.codesigning.azure.net/` |

See [GitHub environment configuration](https://docs.github.com/en/actions/how-tos/deploy/configure-and-manage-deployments/manage-environments)
and the [Azure signing action](https://github.com/Azure/artifact-signing-action).

## Verify a build

Push the workflow changes to an allowed branch, or run **Actions > Build Plugins >
Run workflow** on an allowed branch containing the changes. Inspect the
`sign-windows` job. Its summary records the publisher, certificate thumbprint, and
SHA-256 of the signed plugin. Download `ctools_windows` from that run to test it
in Stata on Windows with application control enabled. Restart Stata before loading
a replaced plugin. Do not install the `windows-unsigned` intermediate artifact.

If Azure login fails, check the exact federation subject, tenant, application ID,
and environment spelling. If signing returns HTTP 403, check that the application
has the signing role on the account containing the active Public Trust profile,
and that the endpoint matches its region. A successful signing job alone does
not imply that the other platform builds or release publication have passed.
