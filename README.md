# QAOA Compiler with Qiskit Backend

This repository includes implementations of QAOA-specific compilation policies described in the following articles (on top of the Qiskit compiler backend):
* [Circuit Compilation Methodologies for Quantum Approximate Optimization Algorithm [MICRO2020]](https://ieeexplore.ieee.org/iel7/9251289/9251849/09251960.pdf?casa_token=bS5P0P7L7W8AAAAA:dWCmbrvC98xyVninXaiZoDweR1C3tPJ7q9YCTeLq_1SPa_HBC6-GBdHjGwgvgigid8CDjoA)
* [An efficient circuit compilation flow for quantum approximate optimization algorithm [DAC2020]](https://ieeexplore.ieee.org/iel7/9211868/9218488/09218558.pdf?casa_token=OkGG7zPyUMYAAAAA:qcKOGwF3jq3tW9F4bfVBoYW78tDlGZnD4LhjXhLN51kreccFfOqiQLEDdvYtRgHQFabAGPI)
* [Noise resilient compilation policies for quantum approximate optimization algorithm [ICCAD2020]](https://dl.acm.org/doi/pdf/10.1145/3400302.3415745?casa_token=D_dOwFq1iIsAAAAA:8mS78EK6GYdV7ELjeh01mi-3lSZRgI9yWeWtYq2o5VBHiCooCPFGZDI5PVbcE12ezLOGNOBDno4)

## Dependencies
* **python>=3.8.5**
* **qiskit>=0.23.1**
* **networkx>=2.5**


## Structure
* [CompileQAOA](): Python Class implementation of the compilation policies.
* [Examples](): Sample input files to the compiler.

## How to Run
```
python main.py {-arg} {arg_val}
```
* {-device_json} {string} : Target device configuration file in json format (basis gates, reliability, allowed two-qubit operations, etc.). An example can be found under [Examples]()
