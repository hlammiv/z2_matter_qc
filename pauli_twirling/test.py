import circuit_twirling.paulitwirling

from qiskit import QuantumCircuit
from qiskit.quantum_info.operators import Operator

if __name__ == "__main__":

    print("testing twirl_cnot" + "." * 5)
    qc = QuantumCircuit(2)
    qc.cx(0, 1)

    check_flag = True
    for i in range(20):
        qc2 = QuantumCircuit(2)
        qc2 = circuit_twirling.paulitwirling.twirl_cnot(qc2, 0, 1)
        check_flag = Operator(qc).equiv(qc2)
        if check_flag is False:
            break
    if check_flag:
        print("twirl_cnot appears to be working")
    else:
        print("twirl_cnot is not working; not equivalent to CNOT")

    print('==' * 50)
    qc = QuantumCircuit(5)
    qc.cx(0, 1)
    qc.cx(2, 3)
    qc2 = QuantumCircuit(5)
    pairs = [(0, 1), (2, 3)]
    if Operator(qc2).equiv(qc):
        print('twirl_hard_cycle is working')
    else:
        print('twirl_hard_cycle is not working')
