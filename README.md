# QAOA Compiler with Qiskit Backend

This repository includes implementations of QAOA-specific compilation policies described in the following articles (on top of the [Qiskit](https://github.com/Qiskit/qiskit) compiler backend):
* [Circuit Compilation Methodologies for Quantum Approximate Optimization Algorithm [MICRO2020]](https://ieeexplore.ieee.org/iel7/9251289/9251849/09251960.pdf?casa_token=bS5P0P7L7W8AAAAA:dWCmbrvC98xyVninXaiZoDweR1C3tPJ7q9YCTeLq_1SPa_HBC6-GBdHjGwgvgigid8CDjoA)
* [An efficient circuit compilation flow for quantum approximate optimization algorithm [DAC2020]](https://ieeexplore.ieee.org/iel7/9211868/9218488/09218558.pdf?casa_token=OkGG7zPyUMYAAAAA:qcKOGwF3jq3tW9F4bfVBoYW78tDlGZnD4LhjXhLN51kreccFfOqiQLEDdvYtRgHQFabAGPI)
* [Noise resilient compilation policies for quantum approximate optimization algorithm [ICCAD2020]](https://dl.acm.org/doi/pdf/10.1145/3400302.3415745?casa_token=D_dOwFq1iIsAAAAA:8mS78EK6GYdV7ELjeh01mi-3lSZRgI9yWeWtYq2o5VBHiCooCPFGZDI5PVbcE12ezLOGNOBDno4)

## Dependencies
* **python>=3.8.5**
* **qiskit>=0.23.1**
* **networkx>=2.5**
* **commentjson>=0.9.0**


## Structure
* [CompileQAOA](https://github.com/mahabubul-alam/QAOA_Compiler/tree/main/CompileQAOA): Python Class implementation of the compilation policies.
* [Examples](https://github.com/mahabubul-alam/QAOA_Compiler/tree/main/Examples): Sample input files to the compiler.

## How to Run
```
python main.py -arg arg_val
```
* -device_json string (mandatory): Target device configuration file location. This file holds the information on basis gates, reliability, and allowed two-qubit operations. It has to be written in json format. An example can be found [here](https://github.com/mahabubul-alam/QAOA_Compiler/blob/main/Examples/QC.json).

* -circuit_json string (mandatory): Problem QAOA-circuit file location. This file holds the required ZZ interactions between various qubit-pairs to encode the cost hamiltonian. It has to be written in json format. An example can be found [here](https://github.com/mahabubul-alam/QAOA_Compiler/blob/main/Examples/QAOA_circ.json).

* -config_json string (mandatory): Compiler configuration file location. This file holds target p-level, and chosen packing limit, qiskit transpiler seed, optimization level, and routing method. It has to be written in json format. An example can be found [here](https://github.com/mahabubul-alam/QAOA_Compiler/blob/main/Examples/Config.json).

* -policy_compilation string: Chosen compilation policy. The current version supports the following policies: Instruction Parallelization-only ('IP'), Iterative Compilation ('IterC'), Incremental Compilation ('IC'), Variation-aware Incremental Compilation ('VIC'). The default value is 'IC'.

* -target_IterC string: Minimization objective of Iterative Compilation. The current version supports the following minimization objectives: Circuit Depth ('D'), Native two-qubit gate-count ('GC_2Q'), Estimated Success Probability ('ESP'). The default value is 'GC_2Q'.

* -output_qasm_file_name string: File name to write the compiled parametric QAOA circuit. The output is written in qasm format. The default value is 'QAOA.qasm'.

### Run the Example
```
python main.py -device_json Examples/QC.json -circuit_json Examples/QAOA_circ.json -config_json Examples/Config.json  -policy_compilation VIC
```
```
python main.py -d Examples/QC.json -ci Examples/QAOA_circ.json -co Examples/Config.json  -p VIC
```

### Notes on the Output File
* The circuit is written with beta/gamma parameters at different P-lavels (https://arxiv.org/pdf/1411.4028.pdf)
* bX --> beta parameter at p=X
* gX_Y_Z --> gamma parameter at p=X between *logical* qubit Y and Z. For the MaxCut problem of unweighted graphs, gX_Y1_Z1 = gX_Y2_Z2 (https://arxiv.org/pdf/1411.4028.pdf)


## Citation
```
@inproceedings{alam2020circuit,
  title={Circuit Compilation Methodologies for Quantum Approximate Optimization Algorithm},
  author={Alam, Mahabubul and Ash-Saki, Abdullah and Ghosh, Swaroop},
  booktitle={2020 53rd Annual IEEE/ACM International Symposium on Microarchitecture (MICRO)},
  pages={215--228},
  year={2020},
  organization={IEEE}
}

@inproceedings{alam2020efficient,
  title={An efficient circuit compilation flow for quantum approximate optimization algorithm},
  author={Alam, Mahabubul and Ash-Saki, Abdullah and Ghosh, Swaroop},
  booktitle={2020 57th ACM/IEEE Design Automation Conference (DAC)},
  pages={1--6},
  year={2020},
  organization={IEEE}
}

@inproceedings{alam2020noise,
  title={Noise resilient compilation policies for quantum approximate optimization algorithm},
  author={Alam, Mahabubul and Ash-Saki, Abdullah and Li, Junde and Chattopadhyay, Anupam and Ghosh, Swaroop},
  booktitle={Proceedings of the 39th International Conference on Computer-Aided Design},
  pages={1--7},
  year={2020}
}
```
