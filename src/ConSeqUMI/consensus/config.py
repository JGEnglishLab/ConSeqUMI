from typing import List

LAST_TRAIN_PATH = {"ltp":"[PATH TO LAST-TRAIN FILE]"}
lamassembleCommandLine = f"lamassemble mat_path_filled_in_programmatically --end -g60 -m 40"

LCOMMAND: List[str] = lamassembleCommandLine.split()

medakaCommandLine = "medaka_consensus -f -m r941_min_high_g303"
MCOMMAND = medakaCommandLine.split()

