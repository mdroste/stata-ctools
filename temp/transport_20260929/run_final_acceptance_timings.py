from pathlib import Path
import subprocess,sys
r=Path(__file__).resolve().parent
for command in ([sys.executable,str(r/'run_campaign.py'),str(r/'final_confirmed')],[sys.executable,str(r/'run_dual.py')],[sys.executable,str(r/'run_commands_final.py')]):
 subprocess.run(command,check=True)
