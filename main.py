import argparse
import os
from CompileQAOA import Compile_QAOA_Qiskit as CQQ


parser = argparse.ArgumentParser(allow_abbrev=True)
parser.add_argument("-device_json", help="Target Device Configuration File", type=str, default=None, action='store', dest='QCD')
parser.add_argument("-circuit_json", help="Problem QAOA ZZ Interaction File", type=str, default=None, action='store', dest='CKT')
parser.add_argument("-config_json", help="QAOA/Compiler Configuration File", type=str, default=None, action='store', dest='CON')
parser.add_argument("-policy_compilation", help="Chosen Compilation Policy: Instruction Parallelization-only (IP), Iterative Compilation (IterC), Incremental Compilation (IC), Variation-aware Incremental Compilation (VIC)", type=str, default='IC', action='store', dest='POL')
parser.add_argument("-target_IterC", help="Iterative Compilation (IterC) should minimize: Depth (D) or Two-qubit gate-count (GC_2Q) or Estimated Success Probability (ESP)", type=str, default='GC_2Q', action='store', dest='TAR')
parser.add_argument("-output_qasm_file_name", help="Name of the output file to write compiled circuit description in QASM format", type=str, default='QAOA.qasm', action='store', dest='OUT')
ARGS = parser.parse_args()




if __name__ == '__main__':
    if os.path.isfile(ARGS.QCD) and os.path.isfile(ARGS.CKT) and os.path.isfile(ARGS.CON):
        Comp_obj = CQQ.Compile_QAOA_Qiskit(Circuit_json=ARGS.CKT, QC_json=ARGS.QCD, Config_json=ARGS.CON, Out_Circuit_File_Name=ARGS.OUT)
        if ARGS.POL == 'IP': Comp_obj.run_IP()
        elif ARGS.POL == 'IterC': Comp_obj.run_IterC(ARGS.TAR)
        elif ARGS.POL == 'IC': Comp_obj.run_IncrC()
        elif ARGS.POL == 'VIC': Comp_obj.run_IncrC(Variation_Aware=True)


