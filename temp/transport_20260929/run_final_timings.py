from pathlib import Path
import subprocess, sys
root = Path(__file__).resolve().parent
subprocess.run([sys.executable, str(root / 'run_campaign.py'), str(root / 'final_revision')], check=True)
subprocess.run([sys.executable, str(root / 'run_commands.py')], check=True)
