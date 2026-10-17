#!/usr/bin/env python3
"""Set up HTTPS access so an AI coding agent can run and inspect Jenkins tests.

The user runs this script in an interactive terminal and supplies their own
Jenkins username and API token, a credential generated in their account settings.
The saved token lets an agent use the Jenkins API (its interface for programs)
without clicking through a browser or logging into the server through SSH.
This script sets up access; the clubb-jenkins skill describes job operations.

Usage:
    From the repository root: python3 jenkins_tests/jenkins_api_login.py
    Generate a named token in your Jenkins account's Security/Configure page,
    then enter your username and token at the prompts. Token input is hidden.

    Public routines:
    - main(argv=None): run interactive setup and return an exit status.
    - verify_credentials(credentials): check a USERNAME:API_TOKEN string over
      HTTPS and return the authenticated username without saving the token.
    - save_credentials(credentials, path): save an already verified string to
      the supplied pathlib.Path, replacing any previous credential file.

Understanding:
    main() verifies the login before writing ~/.config/jenkins/api-credentials.
    The file holds one USERNAME:API_TOKEN line. It is plaintext, with owner-only
    mode 600 permissions on POSIX systems, rather than encrypted storage.
    Rerunning setup replaces this file only after the new login succeeds.

    Verification calls the configured JENKINS_URL's whoAmI/api/json endpoint
    with a 20-second timeout. It rejects redirects and anonymous identities;
    an HTTP 200 response alone does not establish a user login. No persistent
    browser session is created. Later requests authenticate with the saved token
    and remain subject to that user's Jenkins permissions.

References:
    https://www.jenkins.io/doc/book/system-administration/authenticating-scripted-clients/
    https://www.jenkins.io/doc/book/using/remote-access-api/
"""

import argparse
import base64
import getpass
import json
import os
from pathlib import Path
import sys
import tempfile
import urllib.error
import urllib.request
import warnings


JENKINS_URL = "https://carson.math.uwm.edu/jenkins/"
CREDENTIALS_FILE = Path.home() / ".config" / "jenkins" / "api-credentials"


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    """Reject redirects instead of sending the credential to another endpoint."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        """Return None for urllib's request/response arguments to refuse a redirect."""
        return None


def verify_credentials(credentials):
    """Check a USERNAME:API_TOKEN string and return Jenkins's authenticated name.

    Send an identity request without saving the credential or printing the token.
    Raise ValueError for malformed credentials, HTTP errors, invalid JSON, or a
    missing or anonymous identity; connection errors propagate to the caller.
    """
    username, separator, token = credentials.partition(":")
    if not separator or not username or not token or "\n" in credentials or "\r" in credentials:
        raise ValueError("Expected a Jenkins username and API token on one line.")
    authorization = base64.b64encode(credentials.encode()).decode()
    request = urllib.request.Request(
        JENKINS_URL + "whoAmI/api/json",
        headers={"Authorization": "Basic " + authorization},
    )
    opener = urllib.request.build_opener(_NoRedirect())
    try:
        with opener.open(request, timeout=20) as response:
            identity = json.load(response)
    except urllib.error.HTTPError as error:
        if error.code in (401, 403):
            raise ValueError(
                f"Jenkins rejected the login (HTTP {error.code}); check your username and API token."
            ) from None
        raise ValueError(f"Jenkins identity request failed (HTTP {error.code}).") from None
    except (json.JSONDecodeError, UnicodeDecodeError):
        raise ValueError("Jenkins did not return a JSON identity response.") from None
    # Jenkins can return HTTP 200 while still identifying the caller as anonymous.
    if (
        not isinstance(identity, dict)
        or not identity.get("authenticated")
        or not identity.get("name")
        or identity["name"] == "anonymous"
    ):
        raise ValueError("Jenkins did not authenticate a user; credentials were not saved.")
    return identity["name"]


def save_credentials(credentials, path):
    """Write an already verified credential string to the supplied pathlib.Path.

    Create missing parent directories and store one line ending in a newline.
    On POSIX, the new file has mode 600. Replace the destination only after the
    temporary file is written; failed writes/replacement preserve the old file.
    Return None and propagate filesystem errors to the caller.
    """
    path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    # Use the destination's directory so replacement stays on one filesystem.
    descriptor, temporary = tempfile.mkstemp(prefix=".api-credentials-", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w") as output:
            output.write(credentials + "\n")
        os.replace(temporary, path)
    finally:
        # A failed write or replacement must not leave a second credential file.
        Path(temporary).unlink(missing_ok=True)


def main(argv=None):
    """Prompt, verify and save a login; return 0 on success or 1 on failure.

    argv supplies command-line arguments, or None uses the process arguments.
    Help exits before prompting. Report setup failures without printing the
    token, and save to CREDENTIALS_FILE only after verification succeeds.
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog="Generate an API token in your Jenkins account's Security/Configure page. "
        "A verified login replaces ~/.config/jenkins/api-credentials. "
        "Uses only the Python standard library.",
        add_help=False,
    )
    parser.add_argument("-h", "-help", action="help", help="Show this help message and exit.")
    parser.parse_args(argv)
    print("Jenkins:", JENKINS_URL)
    print("Create a named API token in your account's Security/Configure page.")
    try:
        username = input("Jenkins username: ").strip()
        # Stop if getpass cannot hide input instead of falling back to echoed input.
        with warnings.catch_warnings():
            warnings.simplefilter("error", getpass.GetPassWarning)
            token = getpass.getpass("Jenkins API token (hidden): ").strip()
        credentials = username + ":" + token
        # Keep the existing saved login until Jenkins confirms the new identity.
        identity = verify_credentials(credentials)
        save_credentials(credentials, CREDENTIALS_FILE)
    except getpass.GetPassWarning:
        print("Run this command in an interactive terminal to hide the API token.", file=sys.stderr)
        return 1
    except (EOFError, KeyboardInterrupt):
        print("\nLogin cancelled.", file=sys.stderr)
        return 1
    except (ValueError, OSError, urllib.error.URLError) as error:
        print(f"Login failed: {error}", file=sys.stderr)
        return 1
    print("Verified Jenkins identity:", identity)
    print("Saved credentials in:", CREDENTIALS_FILE)
    return 0


if __name__ == "__main__":
    sys.exit(main())
