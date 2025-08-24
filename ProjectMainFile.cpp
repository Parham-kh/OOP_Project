// ----------------- STANDARD LIBRARIES -----------------
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <map>
#include <fstream>
#include <filesystem>
#include <memory>
#include <algorithm>
#include <complex>
#include <regex>
#include <iomanip>
#include <unordered_map>

// ----------------- VENDOR LIBRARIES -----------------
// Graphics & Windowing
#include <GL/glew.h>
#include <GLFW/glfw3.h>

// Dear ImGui
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

// ImPlot for plotting
#include "implot.h"

// Cereal for Save/Load
#include <cereal/archives/json.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/polymorphic.hpp>

// Portable File Dialogs for Open/Save
#include "portable-file-dialogs.h"

// ----------------- PROJECT HEADERS -----------------
//#include "save_load.h"
//#include "recent_files.h"


// ================ Global Application State ================
// Enum to represent the currently selected tool
enum ToolType { NONE, RESISTOR, CAPACITOR, INDUCTOR, VSOURCE, CSOURCE, GROUND, WIRE, VSIN, CSIN, VOLTMETER, DIODE };
ToolType g_currentTool = NONE;


struct PlacedProbe {
    ImVec2 pos;
    std::string node_name;
    double measured_voltage;
};
std::vector<PlacedProbe> g_probes;


// Struct to hold information about a placed graphical element
struct PlacedElement {
    ToolType type;
    ImVec2 pos; // Using ImGui's vector type for position (center of the element)
    std::string name;
    double value;
    bool isVertical = false;
    double amp = 0.0;
    double freq = 0.0;
    double phase = 0.0;
};
std::vector<PlacedElement> g_placedElements;

// A wire is just a pair of points
struct Wire {
    ImVec2 p1, p2;
};
std::vector<Wire> g_wires;

bool g_showComponentPopup = false;
ToolType g_pendingTool = NONE; // The tool waiting for the dialog
char g_componentNameBuffer[64] = "";
char g_componentValueBuffer[64] = "";

bool g_isPlacingWire = false;
ImVec2 g_wireStartPos;

// --- State for the Run Simulation Dialog ---
bool g_showRunPopup = false;
enum RunChoice { RUN_DC, RUN_AC, RUN_VT, RUN_IT };
RunChoice g_runChoice = RUN_DC;

// Buffers for text inputs in the Run dialog
char g_acStartText[64] = "10";
char g_acStopText[64] = "10k";
char g_acPointsPerDecText[64] = "10";
char g_acProbeText[64] = "n1";

char g_trStartText[64] = "0";
char g_trEndText[64] = "1m";
char g_trMaxStepText[64] = "10u";
char g_trProbeText[64] = "n1";

char g_itStartText[64] = "0";
char g_itEndText[64] = "1m";
char g_itMaxStepText[64] = "10u";
char g_itElemText[64] = "R1";


static const float ELEMENT_HALF_LEN = 30.0f;
static const float PIN_STUB = 10.0f;
static const float PIN_OFFSET = ELEMENT_HALF_LEN + PIN_STUB;

bool g_showDCResults = false;
std::map<std::string, double> g_dcVoltages;
std::map<std::string, double> g_dcCurrents;
std::string g_simulationError;

bool g_showDebugInfo = false;

bool g_previewIsVertical = false;

std::string g_componentPopupError;

bool g_showTransientPlot = false;

bool g_showACPlot = false;
std::vector<double> g_ac_freqs;
std::vector<double> g_ac_mags_db;
std::vector<double> g_ac_phases;
std::string g_ac_probe_node;

bool g_showVACPopup = false;
char g_vacNameBuffer[64] = "V2";
char g_vacDcOffsetText[64] = "0";
char g_vacAmpText[64] = "1";
char g_vacFreqText[64] = "1k";

bool g_showIACPopup = false;
char g_iacNameBuffer[64] = "I2";
char g_iacDcOffsetText[64] = "0";
char g_iacAmpText[64] = "1";
char g_iacFreqText[64] = "1k";

// ================ Simulation Engine and Circuit Logic ================
// This is the core simulation engine ported from your original project.

// Custom Exception Classes
class NoGroundException : public std::exception {
public:
    const char* what() const noexcept override { return "Error: No ground node detected in the circuit."; }
};
class InvalidValueException : public std::exception {
public:
    const char* what() const noexcept override { return "Error: Negative or zero value for a component is invalid."; }
};
class SingularMatrixException : public std::exception {
public:
    const char* what() const noexcept override { return "Error: Zero pivot encountered. The system is singular."; }
};

// Wave storage
struct WaveStore {
    std::vector<double> t;
    std::map<std::string, std::vector<double>> V; // node voltages
    std::map<std::string, std::vector<double>> I; // element currents
    void clear(){ t.clear(); V.clear(); I.clear(); }
};
WaveStore g_waves;

// Gaussian Elimination for solving the matrix
std::vector<double> gaussianElimination(std::vector<std::vector<double>> A, std::vector<double> b) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (abs(A[k][i]) > abs(A[maxRow][i])) {
                maxRow = k;
            }
        }
        std::swap(A[i], A[maxRow]);
        std::swap(b[i], b[maxRow]);

        double pivot = A[i][i];
        if (abs(pivot) < 1e-12)
            throw SingularMatrixException();

        for (int j = i; j < n; j++) A[i][j] /= pivot;
        b[i] /= pivot;

        for (int k = i + 1; k < n; k++) {
            double factor = A[k][i];
            for (int j = i; j < n; j++) A[k][j] -= factor * A[i][j];
            b[k] -= factor * b[i];
        }
    }
    std::vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i + 1; j < n; j++) x[i] -= A[i][j] * x[j];
    }
    return x;
}

// Circuit Classes (Nodes/Elements)
class Node {
    std::string name;
public:
    Node(const std::string& nodeName) : name(nodeName) {}
    std::string getName() const { return name; }
};

class Element {
protected:
    std::string name;
    Node* n1;
    Node* n2;
    double value;
public:
    Element(const std::string &elemName, Node* node1, Node* node2, double val)
            : name(elemName), n1(node1), n2(node2), value(val) {}
    virtual ~Element() {}
    virtual std::string getType() const = 0;
    std::string getName() const { return name; }
    double getValue() const { return value; }
    Node* getNode1() const { return n1; }
    Node* getNode2() const { return n2; }
};

class Resistor : public Element {
public:
    Resistor(const std::string &elemName, Node* a, Node* b, double r) : Element(elemName, a, b, r) {
        if (r < 0) throw InvalidValueException();
    }
    std::string getType() const override { return "Resistor"; }
};

class Capacitor : public Element {
public:
    Capacitor(const std::string &elemName, Node* a, Node* b, double c) : Element(elemName, a, b, c) {
        if (c <= 0) throw InvalidValueException();
    }
    std::string getType() const override { return "Capacitor"; }
};

class Inductor : public Element {
public:
    Inductor(const std::string &elemName, Node* a, Node* b, double l) : Element(elemName, a, b, l) {
        if (l <= 0) throw InvalidValueException();
    }
    std::string getType() const override { return "Inductor"; }
};

class VoltageSource : public Element {
public:
    double dcOffset;
    double amplitude = 0;
    double frequency = 0;
    bool isSine = false;

    VoltageSource(const std::string &elemName, Node* n1, Node* n2, double offset, double amp = 0, double freq = 0)
            : Element(elemName, n1, n2, offset),
              dcOffset(offset),
              amplitude(amp),
              frequency(freq) {
        if (amp > 0 && freq > 0) {
            isSine = true;
        }
    }

    std::string getType() const override { return "VoltageSource"; }

    void setValueAtTime(double time) {
        if (isSine) {
            this->value = dcOffset + amplitude * sin(2 * M_PI * frequency * time);
        } else {
            this->value = dcOffset;
        }
    }
};


class CurrentSource : public Element {
public:
    double dcOffset;
    double amplitude = 0;
    double frequency = 0;
    bool isSine = false;

    CurrentSource(const std::string &elemName, Node* a, Node* b, double offset, double amp = 0, double freq = 0)
            : Element(elemName, a, b, offset), dcOffset(offset), amplitude(amp), frequency(freq) {
        if (amp != 0 && freq != 0) {
            isSine = true;
        }
    }

    std::string getType() const override { return "CurrentSource"; }

    void setValueAtTime(double time) {
        if (isSine) {
            this->value = dcOffset + amplitude * sin(2 * M_PI * frequency * time);
        } else {
            this->value = dcOffset;
        }
    }
};


class Diode : public Element {
public:
    // The 'value' from the Element base class will store our vOn
    Diode(const std::string &elemName, Node* anode, Node* cathode, double vOn)
            : Element(elemName, anode, cathode, vOn) {}
    std::string getType() const override { return "Diode"; }
};


// Main Circuit Class
class Circuit {
public:
    std::vector<Node*> nodes;
    std::vector<Element*> elements;
    std::string groundName = "";
    double timeStep = 0.001;
    int totalSteps = 1000;

    ~Circuit() { reset(); }

    void reset() {
        for (auto node : nodes) delete node;
        nodes.clear();
        for (auto elem : elements) delete elem;
        elements.clear();
        groundName = "";
    }

    Node* getOrCreateNode(const std::string& nodeName) {
        for (auto node : nodes)
            if (node->getName() == nodeName) return node;
        Node* newNode = new Node(nodeName);
        nodes.push_back(newNode);
        return newNode;
    }

    void addElement(Element* elem) { elements.push_back(elem); }
    void setGroundNode(const std::string& nodeName) { groundName = nodeName; }

    void setTransientWindow(double t0, double t1, double maxdt) {
        if (t1 <= t0 || maxdt <= 0) throw InvalidValueException();
        int steps = (int)ceil((t1 - t0) / maxdt);
        if (steps < 1) steps = 1;
        this->timeStep = (t1 - t0) / steps;
        this->totalSteps = steps;
    }

    void simulateTransientCapture(WaveStore& ws, const std::vector<std::string>& wanted_nodes, const std::vector<std::string>& wanted_elements) {
        ws.clear();
        double tStart = 0.0;
        double tEnd = (this->totalSteps) * this->timeStep;
        double step = this->timeStep;

        // --- State is now managed locally here, not in the Element classes ---
        std::map<std::string, double> node_voltages;
        std::map<std::string, double> capacitor_prev_voltages;
        std::map<std::string, double> inductor_prev_currents;
        std::map<std::string, double> element_currents;

        // --- Setup ---
        std::map<std::string, int> nodeIndex;
        int idx = 0;
        for (Node* n : nodes) {
            if (n->getName() != this->groundName) {
                nodeIndex[n->getName()] = idx++;
            }
        }

        int n = idx;
        std::vector<VoltageSource*> voltageSources;
        for (Element* e : elements) {
            if (auto vs = dynamic_cast<VoltageSource*>(e)) {
                voltageSources.push_back(vs);
            }
        }
        int m = voltageSources.size();
        int matrixSize = n + m;

        // Initialize stateful components in our local maps
        for (Element* e : elements) {
            if (auto c = dynamic_cast<Capacitor*>(e)) capacitor_prev_voltages[c->getName()] = 0.0;
            else if (auto i = dynamic_cast<Inductor*>(e)) inductor_prev_currents[i->getName()] = 0.0;
        }

        // --- Main Simulation Loop ---
        for (double t = 0; t <= tEnd + 1e-12; t += step) {
            std::vector<std::vector<double>> G(n, std::vector<double>(n, 0.0));
            std::vector<std::vector<double>> B(n, std::vector<double>(m, 0.0));
            std::vector<std::vector<double>> D(m, std::vector<double>(n, 0.0));
            std::vector<double> I(n, 0.0);
            std::vector<double> F(m, 0.0);

            // --- Stamping ---
            for (Element* e : elements) {
                Node* pn = e->getNode1();
                Node* nn = e->getNode2();
                int pi = (pn->getName() != this->groundName) ? nodeIndex.at(pn->getName()) : -1;
                int ni = (nn->getName() != this->groundName) ? nodeIndex.at(nn->getName()) : -1;

                if (auto r = dynamic_cast<Resistor*>(e)) {
                    double g = 1.0 / std::max(r->getValue(), 1e-9);
                    if (pi != -1) G[pi][pi] += g;
                    if (ni != -1) G[ni][ni] += g;
                    if (pi != -1 && ni != -1) { G[pi][ni] -= g; G[ni][pi] -= g; }
                } else if (auto c = dynamic_cast<Capacitor*>(e)) {
                    double geq = c->getValue() / step;
                    double ieq = geq * capacitor_prev_voltages.at(c->getName());
                    if (pi != -1) G[pi][pi] += geq;
                    if (ni != -1) G[ni][ni] += geq;
                    if (pi != -1 && ni != -1) { G[pi][ni] -= geq; G[ni][pi] -= geq; }
                    if (pi != -1) I[pi] += ieq;
                    if (ni != -1) I[ni] -= ieq;
                } else if (auto ind = dynamic_cast<Inductor*>(e)) {
                    double geq = step / std::max(ind->getValue(), 1e-12);
                    double prevI = inductor_prev_currents.at(ind->getName());
                    if (pi != -1) G[pi][pi] += geq;
                    if (ni != -1) G[ni][ni] += geq;
                    if (pi != -1 && ni != -1) { G[pi][ni] -= geq; G[ni][pi] -= geq; }
                    if (pi != -1) I[pi] -= prevI;
                    if (ni != -1) I[ni] += prevI;
                }
                else if (auto cs = dynamic_cast<CurrentSource*>(e)) {
                    cs->setValueAtTime(t); // Update to its value at the current time
                    double I_inst = cs->getValue();
                    if (pi != -1) I[pi] -= I_inst; // Current leaves the positive node
                    if (ni != -1) I[ni] += I_inst; // Current enters the negative node
                }
            }

            for (int k=0; k<m; ++k) {
                auto vs = voltageSources[k];
                Node* pn = vs->getNode1();
                Node* nn = vs->getNode2();
                int pi = (pn->getName() != this->groundName) ? nodeIndex.at(pn->getName()) : -1;
                int ni = (nn->getName() != this->groundName) ? nodeIndex.at(nn->getName()) : -1;

                vs->setValueAtTime(t);
                if (pi != -1) { B[pi][k] = 1.0; D[k][pi] = 1.0; }
                if (ni != -1) { B[ni][k] = -1.0; D[k][ni] = -1.0; }
                F[k] = vs->getValue();
            }

            std::vector<std::vector<double>> A(matrixSize, std::vector<double>(matrixSize, 0.0));
            std::vector<double> Z(matrixSize, 0.0);
            // ... (Assemble A and Z matrices - unchanged)
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) A[i][j] = G[i][j];
                for (int j = 0; j < m; j++) A[i][n + j] = B[i][j];
                Z[i] = I[i];
            }
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) A[n + i][j] = D[i][j];
                Z[n + i] = F[i];
            }

            // --- Solve ---
            std::vector<double> solution = gaussianElimination(A, Z);

            // --- Update and Store Results ---
            node_voltages[this->groundName] = 0.0;
            for (const auto& pair : nodeIndex) {
                node_voltages[pair.first] = solution.at(pair.second);
            }

            ws.t.push_back(t);
            for(const auto& name : wanted_nodes) {
                ws.V[name].push_back(node_voltages.at(name));
            }

            for (int k = 0; k < m; k++) {
                element_currents[voltageSources[k]->getName()] = solution.at(n + k);
            }

            for (Element* e : elements) {
                double v = node_voltages.at(e->getNode1()->getName()) - node_voltages.at(e->getNode2()->getName());
                if (auto c = dynamic_cast<Capacitor*>(e)) {
                    element_currents[c->getName()] = (c->getValue() / step) * (v - capacitor_prev_voltages.at(c->getName()));
                    capacitor_prev_voltages[c->getName()] = v;
                } else if (auto ind = dynamic_cast<Inductor*>(e)) {
                    double current = inductor_prev_currents.at(ind->getName()) + (step / ind->getValue()) * v;
                    element_currents[ind->getName()] = current;
                    inductor_prev_currents[ind->getName()] = current;
                } else if (dynamic_cast<Resistor*>(e)) {
                    element_currents[e->getName()] = v / e->getValue();
                }

                for(const auto& name : wanted_elements) {
                    if(e->getName() == name) {
                        ws.I[name].push_back(element_currents.at(name));
                    }
                }
            }
        }
    }


    bool computeDCOP(std::map<std::string, double>& outV, std::map<std::string, double>& outI) const {
        if (groundName.empty()) throw NoGroundException();

        // --- Part 1: Initial Setup for Iterative Solver ---
        std::vector<Diode*> diodes;
        for (auto elem : elements) {
            if (auto d = dynamic_cast<Diode*>(elem)) diodes.push_back(d);
        }
        std::map<std::string, bool> diode_is_on;
        for (auto d : diodes) diode_is_on[d->getName()] = false;

        outV.clear(); outI.clear();

        // --- FIX: Variables to store the final solution state ---
        std::vector<double> final_x;
        int final_n = 0, final_m = 0;
        std::vector<const Element*> final_v_sources_list;

        // --- Part 2: Main Iteration Loop ---
        const int MAX_ITERATIONS = 100;
        for (int iter = 0; iter < MAX_ITERATIONS; ++iter) {
            std::vector<const Element*> v_sources_and_on_diodes;
            for (auto elem : elements) {
                if (dynamic_cast<const VoltageSource*>(elem)) {
                    v_sources_and_on_diodes.push_back(elem);
                } else if (auto d = dynamic_cast<const Diode*>(elem)) {
                    if (diode_is_on.at(d->getName())) {
                        v_sources_and_on_diodes.push_back(d);
                    }
                }
            }

            std::vector<Node*> unknown_nodes;
            for (auto node : nodes) if (node->getName() != groundName) unknown_nodes.push_back(node);

            const int n = unknown_nodes.size();
            const int m = v_sources_and_on_diodes.size();
            const int matrix_size = n + m;
            if (matrix_size == 0) return true;

            std::map<std::string, int> node_idx;
            for (int i = 0; i < n; ++i) node_idx[unknown_nodes[i]->getName()] = i;

            std::vector<std::vector<double>> A(matrix_size, std::vector<double>(matrix_size, 0.0));
            std::vector<double> b(matrix_size, 0.0);

            // --- C. Stamp Linear Elements and OFF Diodes ---
            for (auto elem : elements) {
                std::string name_a = elem->getNode1()->getName();
                std::string name_c = elem->getNode2()->getName();
                bool is_a_unknown = (name_a != groundName);
                bool is_c_unknown = (name_c != groundName);
                int idx_a = is_a_unknown ? node_idx.at(name_a) : -1;
                int idx_c = is_c_unknown ? node_idx.at(name_c) : -1;

                double g = 0;
                if (dynamic_cast<Resistor*>(elem)) { g = 1.0 / std::max(elem->getValue(), 1e-9); }
                else if (dynamic_cast<Capacitor*>(elem)) { g = 1e-12; }
                else if (dynamic_cast<Inductor*>(elem)) { g = 1e9; }
                else if (auto* d = dynamic_cast<Diode*>(elem)) {
                    // --- FIX: Handle OFF diodes ---
                    if (!diode_is_on.at(d->getName())) { g = 1e-12; } // Treat as open circuit
                }
                else if (auto* cs = dynamic_cast<CurrentSource*>(elem)) {
                    double I = cs->getValue();
                    if (is_a_unknown) b[idx_a] -= I;
                    if (is_c_unknown) b[idx_c] += I;
                }

                if (g != 0) {
                    if (is_a_unknown) A[idx_a][idx_a] += g;
                    if (is_c_unknown) A[idx_c][idx_c] += g;
                    if (is_a_unknown && is_c_unknown) {
                        A[idx_a][idx_c] -= g;
                        A[idx_c][idx_a] -= g;
                    }
                }
            }

            // --- D. Stamp Voltage Sources and ON Diodes using MNA ---
            for (int k = 0; k < m; ++k) {
                const auto* elem = v_sources_and_on_diodes[k];
                std::string name_a = elem->getNode1()->getName(); // Positive terminal
                std::string name_c = elem->getNode2()->getName();
                double V = elem->getValue(); // Works for both VSource.value and Diode.vOn
                int current_idx = n + k;

                if (name_a != groundName) { A[node_idx.at(name_a)][current_idx] += 1.0; A[current_idx][node_idx.at(name_a)] += 1.0; }
                if (name_c != groundName) { A[node_idx.at(name_c)][current_idx] -= 1.0; A[current_idx][node_idx.at(name_c)] -= 1.0; }
                b[current_idx] += V;
            }

            std::vector<double> x = gaussianElimination(A, b);
            for (auto& kv : node_idx) outV[kv.first] = x[kv.second];

            // --- F. Check for Convergence ---
            bool converged = true;
            for (auto d : diodes) {
                double Va = outV.count(d->getNode1()->getName()) ? outV.at(d->getNode1()->getName()) : 0.0;
                double Vc = outV.count(d->getNode2()->getName()) ? outV.at(d->getNode2()->getName()) : 0.0;
                bool should_be_on = (Va - Vc >= d->getValue());

                if (should_be_on != diode_is_on.at(d->getName())) {
                    converged = false;
                    diode_is_on[d->getName()] = should_be_on;
                }
            }

            // --- FIX: Save the state of the final successful iteration ---
            final_x = x;
            final_n = n;
            final_m = m;
            final_v_sources_list = v_sources_and_on_diodes;

            if (converged && iter > 0) break;
            if (iter == MAX_ITERATIONS - 1) throw std::runtime_error("DC solution did not converge.");
        }

        // --- Part 3: Final Current Calculation ---
        // --- FIX: Use the saved final_... variables ---
        for (int k = 0; k < final_m; ++k) {
            outI[final_v_sources_list[k]->getName()] = final_x[final_n + k];
        }

        for (auto elem : elements) {
            if (outI.count(elem->getName())) continue;

            double Va = (elem->getNode1()->getName() == groundName) ? 0.0 : outV.at(elem->getNode1()->getName());
            double Vc = (elem->getNode2()->getName() == groundName) ? 0.0 : outV.at(elem->getNode2()->getName());

            if (dynamic_cast<Resistor*>(elem)) {
                outI[elem->getName()] = (Va - Vc) / std::max(elem->getValue(), 1e-9);
            } else if (dynamic_cast<Inductor*>(elem)) {
                // In DC, I = V / R_short = V * G_short
                outI[elem->getName()] = (Va - Vc) * 1e9;
            } else if (dynamic_cast<Capacitor*>(elem)) {
                outI[elem->getName()] = 0.0;
            } else if (dynamic_cast<CurrentSource*>(elem)) {
                outI[elem->getName()] = elem->getValue();
            } else if (dynamic_cast<Diode*>(elem)) {
                // This case handles diodes that were OFF in the final solution.
                // Their current is 0. ON diodes were handled above.
                outI[elem->getName()] = 0.0;
            } else {
                outI[elem->getName()] = NAN;
            }
        }
        return true;
    }


};


Circuit g_circuit; // A global circuit object for the simulation


std::vector<double> logspace(double start, double stop, int num_points) {
    if (start <= 0 || stop <= 0 || stop < start || num_points < 2) {
        return {};
    }
    std::vector<double> points;
    points.reserve(num_points);
    double log_start = log10(start);
    double log_stop = log10(stop);
    double log_step = (log_stop - log_start) / (num_points - 1);
    for (int i = 0; i < num_points; ++i) {
        points.push_back(pow(10.0, log_start + i * log_step));
    }
    return points;
}

// The main AC solver function
using cd = std::complex<double>;
static const cd j(0.0, 1.0);

bool solveACAtFrequency(const Circuit& C, double freq, std::map<std::string, cd>& outV) {
    outV.clear();
    if (freq <= 0.0 || C.groundName.empty()) return false;

    std::vector<Node*> unknown_nodes;
    for (auto* n : C.nodes) if (n->getName() != C.groundName) unknown_nodes.push_back(n);

    std::vector<const VoltageSource*> v_sources;
    for (auto* e : C.elements) if (auto* v = dynamic_cast<const VoltageSource*>(e)) v_sources.push_back(v);

    const int n = unknown_nodes.size();
    const int m = v_sources.size();
    const int matrix_size = n + m;
    if (matrix_size == 0) return true;

    std::map<std::string, int> node_idx;
    for (int i = 0; i < n; ++i) node_idx[unknown_nodes[i]->getName()] = i;

    const double w = 2.0 * M_PI * freq;
    std::vector<std::vector<cd>> A(matrix_size, std::vector<cd>(matrix_size, {0,0}));
    std::vector<cd> b(matrix_size, {0,0});

    // Stamp elements with complex impedance
    for (auto* e : C.elements) {
        std::string name_a = e->getNode1()->getName();
        std::string name_c = e->getNode2()->getName();
        int idx_a = (name_a != C.groundName) ? node_idx.at(name_a) : -1;
        int idx_c = (name_c != C.groundName) ? node_idx.at(name_c) : -1;

        cd y = {0,0}; // Admittance
        if (auto r = dynamic_cast<Resistor*>(e)) {
            y = 1.0 / std::max(r->getValue(), 1e-9);
        }
        else if (auto c = dynamic_cast<Capacitor*>(e)) {
            y = j * w * c->getValue();
        }
        else if (auto l = dynamic_cast<Inductor*>(e)) {
            y = 1.0 / (j * w * std::max(l->getValue(), 1e-12));
        }
        else if (auto cs = dynamic_cast<CurrentSource*>(e)) {
            // For AC, stamp the source's amplitude into the RHS 'b' vector
            double I_amp = cs->amplitude;
            if (idx_a != -1) b[idx_a] -= I_amp; // Current leaves the positive node
            if (idx_c != -1) b[idx_c] += I_amp; // Current enters the negative node
        }

        // Stamp the admittance (if any) into the A matrix
        if (y.real() != 0 || y.imag() != 0) {
            if (idx_a != -1) A[idx_a][idx_a] += y;
            if (idx_c != -1) A[idx_c][idx_c] += y;
            if (idx_a != -1 && idx_c != -1) { A[idx_a][idx_c] -= y; A[idx_c][idx_a] -= y; }
        }
    }

    // Stamp voltage sources for MNA
    for (int k = 0; k < m; ++k) {
        auto vs = v_sources[k];
        int idx_a = (vs->getNode1()->getName() != C.groundName) ? node_idx.at(vs->getNode1()->getName()) : -1;
        int idx_c = (vs->getNode2()->getName() != C.groundName) ? node_idx.at(vs->getNode2()->getName()) : -1;
        int current_idx = n + k;

        if (idx_a != -1) { A[idx_a][current_idx] += 1.0; A[current_idx][idx_a] += 1.0; }
        if (idx_c != -1) { A[idx_c][current_idx] -= 1.0; A[current_idx][idx_c] -= 1.0; }

        // For AC, we use the sine amplitude, not the DC offset
        b[current_idx] = vs->amplitude;
    }

    // Simple complex solver using Gaussian elimination
    auto solveComplex = [&](std::vector<std::vector<cd>> M, std::vector<cd> v)->std::vector<cd> {
        const int size = M.size();
        for (int i = 0; i < size; ++i) {
            int maxRow = i;
            for (int k = i + 1; k < size; ++k) if (std::abs(M[k][i]) > std::abs(M[maxRow][i])) maxRow = k;
            std::swap(M[i], M[maxRow]); std::swap(v[i], v[maxRow]);
            cd p = M[i][i];
            if (std::abs(p) < 1e-12) throw SingularMatrixException();
            for (int k = i; k < size; ++k) M[i][k] /= p;
            v[i] /= p;
            for (int k = 0; k < size; ++k) {
                if (k != i) {
                    cd f = M[k][i];
                    for (int j = i; j < size; ++j) M[k][j] -= f * M[i][j];
                    v[k] -= f * v[i];
                }
            }
        }
        return v;
    };

    std::vector<cd> x = solveComplex(A, b);
    for (int i = 0; i < n; ++i) outV[unknown_nodes[i]->getName()] = x[i];
    return true;
}


// --- Logic to convert graphical elements into a simulatable circuit ---
struct Connector { int id; ImVec2 pos; };
std::vector<Connector> g_allConnectors;
std::map<int, std::vector<int>> g_connectorGraph;
std::map<int, std::string> g_connectorToNodeName;
int g_connectorIdCounter = 0;

int addConnector(ImVec2 pos) {
    for (const auto& c : g_allConnectors) {
        float dx = c.pos.x - pos.x;
        float dy = c.pos.y - pos.y;
        if (dx * dx + dy * dy < 4.0f) return c.id;
    }
    g_allConnectors.push_back({g_connectorIdCounter, pos});
    return g_connectorIdCounter++;
}


void buildCircuit() {
    // 1. Reset everything for a fresh build
    g_circuit.reset();
    g_allConnectors.clear();
    g_connectorGraph.clear();
    g_connectorToNodeName.clear();
    g_connectorIdCounter = 0;

    for (const auto& e : g_placedElements) {
        if (e.type == GROUND) {
            addConnector(e.pos);
        } else {
            if (e.isVertical) {
                addConnector(ImVec2(e.pos.x, e.pos.y - PIN_OFFSET));
                addConnector(ImVec2(e.pos.x, e.pos.y + PIN_OFFSET));
            } else {
                addConnector(ImVec2(e.pos.x - PIN_OFFSET, e.pos.y));
                addConnector(ImVec2(e.pos.x + PIN_OFFSET, e.pos.y));
            }
        }
    }

    for (const auto& w : g_wires) {
        int id1 = addConnector(w.p1);
        int id2 = addConnector(w.p2);
        g_connectorGraph[id1].push_back(id2);
        g_connectorGraph[id2].push_back(id1);
    }

    // 3. Find connected groups (nodes) using Disjoint Set Union
    std::map<int, int> parent;
    std::function<int(int)> find =
            [&](int i) -> int {
                if (parent.find(i) == parent.end() || parent[i] == i) return i;
                return parent[i] = find(parent[i]);
            };
    auto unite = [&](int i, int j) {
        int root_i = find(i);
        int root_j = find(j);
        if (root_i != root_j) parent[root_i] = root_j;
    };

    for (const auto& c : g_allConnectors) parent[c.id] = c.id;
    for (const auto& pair : g_connectorGraph) {
        for (int neighbor : pair.second) {
            unite(pair.first, neighbor);
        }
    }

    // 4. Find the ground node's group ID FIRST
    int ground_root_id = -1;
    for (const auto& e : g_placedElements) {
        if (e.type == GROUND) {
            int conn_id = addConnector(e.pos);
            ground_root_id = find(conn_id);
            break;
        }
    }

    if (ground_root_id == -1) {
        throw NoGroundException();
    }

    // 5. Assign names to each node group, giving ground a special name "0"
    std::map<int, std::string> root_to_name;
    int node_idx_counter = 1;
    for (const auto& c : g_allConnectors) {
        int root = find(c.id);
        if (root_to_name.find(root) == root_to_name.end()) {
            if (root == ground_root_id) {
                root_to_name[root] = "0"; // Special name for ground
            } else {
                root_to_name[root] = "n" + std::to_string(node_idx_counter++);
            }
        }
        g_connectorToNodeName[c.id] = root_to_name[root];
    }

    // 6. Build the logical circuit object
    g_circuit.setGroundNode("0");
    for (const auto& pair : root_to_name) {
        g_circuit.getOrCreateNode(pair.second);
    }

    for (const auto& e : g_placedElements) {
        if (e.type == GROUND || e.type == WIRE || e.type == VOLTMETER) continue;

        int id1, id2;
        if (e.isVertical) {
            id1 = addConnector(ImVec2(e.pos.x, e.pos.y - PIN_OFFSET));
            id2 = addConnector(ImVec2(e.pos.x, e.pos.y + PIN_OFFSET));
        } else {
            id1 = addConnector(ImVec2(e.pos.x - PIN_OFFSET, e.pos.y));
            id2 = addConnector(ImVec2(e.pos.x + PIN_OFFSET, e.pos.y));
        }

        Node* n1 = g_circuit.getOrCreateNode(g_connectorToNodeName.at(id1));
        Node* n2 = g_circuit.getOrCreateNode(g_connectorToNodeName.at(id2));

        switch (e.type) {
            case RESISTOR: g_circuit.addElement(new Resistor(e.name, n1, n2, e.value)); break;
            case CAPACITOR: g_circuit.addElement(new Capacitor(e.name, n1, n2, e.value)); break;
            case INDUCTOR: g_circuit.addElement(new Inductor(e.name, n1, n2, e.value)); break;
            case VSOURCE: g_circuit.addElement(new VoltageSource(e.name, n1, n2, e.value)); break;
            case CSOURCE: g_circuit.addElement(new CurrentSource(e.name, n1, n2, e.value)); break;
            case DIODE: g_circuit.addElement(new Diode(e.name, n1, n2, e.value)); break;
            case VSIN: g_circuit.addElement(new VoltageSource(e.name, n1, n2, e.value, e.amp, e.freq)); break;
            case CSIN: g_circuit.addElement(new CurrentSource(e.name, n1, n2, e.value, e.amp, e.freq)); break;

            default: break;
        }
    }

    // Add wires as ideal resistors
    for (const auto& w : g_wires) {
        int id1 = addConnector(w.p1);
        int id2 = addConnector(w.p2);
        std::string name1 = g_connectorToNodeName.at(id1);
        std::string name2 = g_connectorToNodeName.at(id2);
        if (name1 != name2) {
            Node* n1 = g_circuit.getOrCreateNode(name1);
            Node* n2 = g_circuit.getOrCreateNode(name2);
            g_circuit.addElement(new Resistor("W", n1, n2, 1e-9));
        }
    }
}


// ================ (Add after Global State Variables) ================
double parseNumber(std::string input) {
    if (input.empty()) return 0.0;
    double result = 0;
    char last = input.back();
    std::string numberPart = input;

    if (isalpha(last)) {
        numberPart = input.substr(0, input.length() - 1);
    }

    try {
        result = std::stod(numberPart);
        if (last == 'M') result *= 1e6;
        else if (last == 'k' || last == 'K') result *= 1e3;
        else if (last == 'm') result *= 1e-3;
        else if (last == 'u' || last == 'U') result *= 1e-6;
        else if (last == 'n') result *= 1e-9;
    } catch (...) { return 0.0; }

    return result;
}

std::vector<std::string> splitList(const std::string& s) {
    std::vector<std::string> out;
    std::string tmp;
    std::stringstream ss(s);
    while (ss >> tmp) {
        // Handle comma at the end of a token, e.g., "R1,"
        if (!tmp.empty() && tmp.back() == ',') {
            tmp.pop_back();
        }
        if (!tmp.empty()) {
            out.push_back(tmp);
        }
    }
    return out;
}

// ================ UI Rendering Functions ================
void RenderMenuBar() {
    if (ImGui::BeginMainMenuBar()) {
        if (ImGui::BeginMenu("File")) {
            // ... (File menu items remain the same) ...
            if (ImGui::MenuItem("Exit")) { exit(0); }
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Run")) {
            if (ImGui::MenuItem("Run Simulation...")) {
                g_showRunPopup = true;
            }
            ImGui::EndMenu();
        }
        // --- NEW DEBUG MENU ---
        if (ImGui::BeginMenu("Debug")) {
            ImGui::MenuItem("Show Node Info", NULL, &g_showDebugInfo);
            ImGui::EndMenu();
        }
        ImGui::EndMainMenuBar();
    }
}


void RenderToolbar() {
    ImGui::BeginChild("Toolbar", ImVec2(0, 35), true, ImGuiWindowFlags_NoScrollbar);

    auto request_component_popup = [](ToolType type) {
        g_pendingTool = type;
        g_showComponentPopup = true;
        g_componentNameBuffer[0] = '\0';
        g_componentValueBuffer[0] = '\0';
        g_componentPopupError.clear();
    };

    if (ImGui::Button("R")) { request_component_popup(RESISTOR); } ImGui::SameLine();
    if (ImGui::Button("C")) { request_component_popup(CAPACITOR); } ImGui::SameLine();
    if (ImGui::Button("L")) { request_component_popup(INDUCTOR); } ImGui::SameLine();
    if (ImGui::Button("V")) { request_component_popup(VSOURCE); } ImGui::SameLine();
    if (ImGui::Button("I")) { request_component_popup(CSOURCE); } ImGui::SameLine();
    if (ImGui::Button("D")) { request_component_popup(DIODE); } ImGui::SameLine();
    if (ImGui::Button("GND")) { g_currentTool = GROUND; } ImGui::SameLine();
    if (ImGui::Button("VAC")) { g_showVACPopup = true; } ImGui::SameLine();
    if (ImGui::Button("IAC")) { g_showIACPopup = true; } ImGui::SameLine();
    if (ImGui::Button("VM")) { g_currentTool = VOLTMETER; } ImGui::SameLine();

    if (ImGui::Button("Wire")) {
        g_currentTool = WIRE;
        g_isPlacingWire = false; // This resets the state
    }
    ImGui::SameLine();



    ImGui::EndChild();
}

void RenderComponentPopup() {
    if (g_showComponentPopup) {
        ImGui::OpenPopup("Component Properties");
        g_showComponentPopup = false;
    }

    ImVec2 center = ImGui::GetMainViewport()->GetCenter();
    ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));

    if (ImGui::BeginPopupModal("Component Properties", NULL, ImGuiWindowFlags_AlwaysAutoResize)) {
        ImGui::Text("Enter Component Details");
        ImGui::Separator();

        ImGui::InputText("Name", g_componentNameBuffer, 64);
        ImGui::InputText("Value", g_componentValueBuffer, 64);

        // --- Display error message if it exists ---
        if (!g_componentPopupError.empty()) {
            ImGui::TextColored(ImVec4(1.0f, 0.0f, 0.0f, 1.0f), "%s", g_componentPopupError.c_str());
        }

        if (ImGui::Button("OK", ImVec2(120, 0))) {
            std::string new_name = g_componentNameBuffer;
            bool name_exists = false;

            // --- Validation Logic ---
            if (new_name.empty()) {
                g_componentPopupError = "Error: Name cannot be empty.";
            } else {
                for (const auto& elem : g_placedElements) {
                    if (elem.name == new_name) {
                        name_exists = true;
                        break;
                    }
                }

                if (name_exists) {
                    g_componentPopupError = "Error: Name '" + new_name + "' is already in use.";
                } else {
                    // Success! Arm the tool and close the popup.
                    g_componentPopupError.clear();
                    g_currentTool = g_pendingTool;
                    ImGui::CloseCurrentPopup();
                }
            }
        }
        ImGui::SameLine();
        if (ImGui::Button("Cancel", ImVec2(120, 0))) {
            g_pendingTool = NONE;
            g_currentTool = NONE;
            ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
    }
}


void RenderRunPopup() {
    if (g_showRunPopup) {
        ImGui::OpenPopup("Run Simulation");
        g_showRunPopup = false; // Reset the trigger
    }

    ImVec2 center = ImGui::GetMainViewport()->GetCenter();
    ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));

    if (ImGui::BeginPopupModal("Run Simulation", NULL, ImGuiWindowFlags_AlwaysAutoResize)) {
        ImGui::Text("Select analysis type:");
        ImGui::Separator();

        // --- Radio Buttons for analysis choice ---
        ImGui::RadioButton("DC Operating Point", (int*)&g_runChoice, RUN_DC);
        ImGui::RadioButton("AC Sweep", (int*)&g_runChoice, RUN_AC);
        ImGui::RadioButton("Transient (V vs t)", (int*)&g_runChoice, RUN_VT);
        ImGui::RadioButton("Transient (I vs t)", (int*)&g_runChoice, RUN_IT);
        ImGui::Separator();

        // --- Show inputs based on the selected analysis type ---
        if (g_runChoice == RUN_AC) {
            ImGui::InputText("Start Freq (Hz)", g_acStartText, 64);
            ImGui::InputText("Stop Freq (Hz)", g_acStopText, 64);
            ImGui::InputText("Points/decade", g_acPointsPerDecText, 64);
            ImGui::InputText("Probe Node", g_acProbeText, 64);
        } else if (g_runChoice == RUN_VT) {
            ImGui::InputText("Start Time (s)", g_trStartText, 64);
            ImGui::InputText("End Time (s)", g_trEndText, 64);
            ImGui::InputText("Max Timestep (s)", g_trMaxStepText, 64);
            ImGui::InputText("Probe Node(s)", g_trProbeText, 64);
        } else if (g_runChoice == RUN_IT) {
            ImGui::InputText("Start Time (s)", g_itStartText, 64);
            ImGui::InputText("End Time (s)", g_itEndText, 64);
            ImGui::InputText("Max Timestep (s)", g_itMaxStepText, 64);
            ImGui::InputText("Element(s)", g_itElemText, 64);
        } else {
            ImGui::Text("Calculates the DC operating point of the circuit.");
        }

        ImGui::Separator();
        if (ImGui::Button("Run", ImVec2(120, 0))) {
            g_simulationError.clear();
            try {
                buildCircuit();
                if (g_runChoice == RUN_DC) {
                    g_circuit.computeDCOP(g_dcVoltages, g_dcCurrents);
                    g_showDCResults = true;
                } else if (g_runChoice == RUN_AC) {
                    double f_start = parseNumber(g_acStartText);
                    double f_stop = parseNumber(g_acStopText);
                    int points = std::stoi(g_acPointsPerDecText) * log10(f_stop/f_start); // Pts/decade * num_decades
                    g_ac_probe_node = g_acProbeText;

                    g_ac_freqs = logspace(f_start, f_stop, std::max(10, points));
                    g_ac_mags_db.clear();
                    g_ac_phases.clear();

                    for (double freq : g_ac_freqs) {
                        std::map<std::string, cd> V;
                        solveACAtFrequency(g_circuit, freq, V);

                        cd v_out = V.count(g_ac_probe_node) ? V.at(g_ac_probe_node) : 0;
                        double mag = std::abs(v_out);
                        double mag_db = 20 * log10(std::max(mag, 1e-9)); // Convert to dB
                        double phase = std::arg(v_out) * 180.0 / M_PI;

                        g_ac_mags_db.push_back(mag_db);
                        g_ac_phases.push_back(phase);
                    }
                    g_showACPlot = true;
                } else if (g_runChoice == RUN_VT) {
                    double t_start = 0.0;
                    double t_end = parseNumber(g_trEndText);
                    double t_step = parseNumber(g_trMaxStepText);

                    // Parse the node list and call the function with all arguments
                    std::vector<std::string> wanted_nodes = splitList(g_trProbeText);
                    std::vector<std::string> wanted_elements; // Empty for V-t analysis

                    g_circuit.setTransientWindow(t_start, t_end, t_step);
                    g_circuit.simulateTransientCapture(g_waves, wanted_nodes, wanted_elements);
                    g_showTransientPlot = true;
                } else if (g_runChoice == RUN_IT) {
                    double t_start = 0.0;
                    double t_end = parseNumber(g_itEndText);
                    double t_step = parseNumber(g_itMaxStepText);

                    // Parse the element list and call the function with all arguments
                    std::vector<std::string> wanted_nodes; // Empty for I-t analysis
                    std::vector<std::string> wanted_elements = splitList(g_itElemText);

                    g_circuit.setTransientWindow(t_start, t_end, t_step);
                    g_circuit.simulateTransientCapture(g_waves, wanted_nodes, wanted_elements);
                    g_showTransientPlot = true;
                }
            } catch (const std::exception& e) {
                g_simulationError = e.what();
                g_showDCResults = true;
            }
            ImGui::CloseCurrentPopup();
        }
        ImGui::SameLine();
        if (ImGui::Button("Cancel", ImVec2(120, 0))) {
            ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
    }
}


void RenderDCResultsWindow() {
    if (!g_showDCResults) return;

    ImGui::Begin("DC Operating Point", &g_showDCResults, ImGuiWindowFlags_AlwaysAutoResize);

    // Display any errors first
    if (!g_simulationError.empty()) {
        ImGui::TextColored(ImVec4(1, 0, 0, 1), "Error:");
        ImGui::TextWrapped("%s", g_simulationError.c_str());
    } else {
        // Display Voltages
        ImGui::Text("Nodal Voltages");
        if (ImGui::BeginTable("Voltages", 2, ImGuiTableFlags_Borders)) {
            ImGui::TableSetupColumn("Node");
            ImGui::TableSetupColumn("Voltage (V)");
            ImGui::TableHeadersRow();
            for (const auto& pair : g_dcVoltages) {
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0);
                ImGui::Text("%s", pair.first.c_str());
                ImGui::TableSetColumnIndex(1);
                ImGui::Text("%.4f", pair.second);
            }
            ImGui::EndTable();
        }

        ImGui::Spacing();

        // Display Currents
        ImGui::Text("Element Currents");
        if (ImGui::BeginTable("Currents", 2, ImGuiTableFlags_Borders)) {
            ImGui::TableSetupColumn("Element");
            ImGui::TableSetupColumn("Current (A)");
            ImGui::TableHeadersRow();
            for (const auto& pair : g_dcCurrents) {
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0);
                ImGui::Text("%s", pair.first.c_str());
                ImGui::TableSetColumnIndex(1);
                if (std::isnan(pair.second)) {
                    ImGui::Text("N/A");
                } else {
                    ImGui::Text("%.4f", pair.second);
                }
            }
            ImGui::EndTable();
        }
    }

    if (ImGui::Button("Close")) {
        g_showDCResults = false;
    }
    ImGui::End();
}


void RenderTransientPlotWindow() {
    if (!g_showTransientPlot) return;

    ImGui::SetNextWindowSize(ImVec2(700, 500), ImGuiCond_FirstUseEver);
    ImGui::Begin("Transient Analysis Results", &g_showTransientPlot);

    if (ImPlot::BeginPlot("Waveforms")) {
        ImPlot::SetupAxes("Time (s)", "Value");

        // Plot all requested node voltages
        for (const auto& pair : g_waves.V) {
            ImPlot::PlotLine(pair.first.c_str(), g_waves.t.data(), pair.second.data(), g_waves.t.size());
        }

        // Plot all requested element currents
        for (const auto& pair : g_waves.I) {
            ImPlot::PlotLine(pair.first.c_str(), g_waves.t.data(), pair.second.data(), g_waves.t.size());
        }

        ImPlot::EndPlot();
    }
    ImGui::End();
}


void RenderACPlotWindow() {
    if (!g_showACPlot) return;

    ImGui::SetNextWindowSize(ImVec2(700, 500), ImGuiCond_FirstUseEver);
    ImGui::Begin("AC Sweep Results", &g_showACPlot);

    // Begin the plot with only a title
    if (ImPlot::BeginPlot(("V(" + g_ac_probe_node + ")").c_str())) {
        ImPlot::SetupAxes("Frequency (Hz)", "Magnitude (dB)");

        ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);

        if (!g_ac_freqs.empty()) {
            ImPlot::PlotLine("Magnitude", g_ac_freqs.data(), g_ac_mags_db.data(), g_ac_freqs.size());
        }

        ImPlot::EndPlot();
    }
    ImGui::End();
}


void RenderVACPopup() {
    if (g_showVACPopup) {
        ImGui::OpenPopup("AC Voltage Source");
        g_showVACPopup = false;
    }

    ImVec2 center = ImGui::GetMainViewport()->GetCenter();
    ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));

    if (ImGui::BeginPopupModal("AC Voltage Source", NULL, ImGuiWindowFlags_AlwaysAutoResize)) {
        ImGui::Text("AC Source Properties");
        ImGui::Separator();

        ImGui::InputText("Name", g_vacNameBuffer, 64);
        ImGui::InputText("DC Offset (V)", g_vacDcOffsetText, 64);
        ImGui::InputText("Amplitude (V)", g_vacAmpText, 64);
        ImGui::InputText("Frequency (Hz)", g_vacFreqText, 64);

        // You can add duplicate name validation here just like in RenderComponentPopup

        if (ImGui::Button("OK", ImVec2(120, 0))) {
            g_currentTool = VSIN; // Arm the VSIN tool for placement
            ImGui::CloseCurrentPopup();
        }
        ImGui::SameLine();
        if (ImGui::Button("Cancel", ImVec2(120, 0))) {
            g_currentTool = NONE;
            ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
    }
}


void RenderIACPopup() {
    if (g_showIACPopup) {
        ImGui::OpenPopup("AC Current Source");
        g_showIACPopup = false;
    }

    // ... (This function is identical to RenderVACPopup, just with different labels and titles)
    // You can copy RenderVACPopup and change the title, labels, and buffers to g_iac...
    if (ImGui::BeginPopupModal("AC Current Source", NULL, ImGuiWindowFlags_AlwaysAutoResize)) {
        ImGui::Text("AC Source Properties");
        ImGui::Separator();
        ImGui::InputText("Name", g_iacNameBuffer, 64);
        ImGui::InputText("DC Offset (A)", g_iacDcOffsetText, 64);
        ImGui::InputText("Amplitude (A)", g_iacAmpText, 64);
        ImGui::InputText("Frequency (Hz)", g_iacFreqText, 64);

        if (ImGui::Button("OK", ImVec2(120, 0))) {
            g_currentTool = CSIN; // Arm the CSIN tool
            ImGui::CloseCurrentPopup();
        }
        ImGui::SameLine();
        if (ImGui::Button("Cancel", ImVec2(120, 0))) { ImGui::CloseCurrentPopup(); }
        ImGui::EndPopup();
    }
}


// ========================================================================================================
// Helper function to format a number with SI prefixes (k, M, m, u, n)
std::string formatValueWithSI(double value) {
    if (value == 0.0) return "0";

    const char* prefixes[] = {"p", "n", "u", "m", "", "k", "M", "G"};
    int index = 4; // Index for the empty prefix ""
    double display_val = value;

    if (std::abs(value) >= 1.0) {
        while (index < 7 && std::abs(display_val) >= 1000.0) {
            display_val /= 1000.0;
            index++;
        }
    } else {
        while (index > 0 && std::abs(display_val) < 1.0) {
            display_val *= 1000.0;
            index--;
        }
    }

    // Use ostringstream for clean formatting without trailing zeros
    std::ostringstream oss;
    oss << display_val << prefixes[index];
    return oss.str();
}


// --- Drawing constants ---
void DrawComponent(ImDrawList* drawList, const PlacedElement& element, const ImVec2& canvas_p0) {
    const ImU32 color = IM_COL32(0, 0, 0, 255);
    const float thickness = 1.5f;
    const ImVec2 pos = ImVec2(canvas_p0.x + element.pos.x, canvas_p0.y + element.pos.y);

    // Calculate start/end points based on rotation
    ImVec2 p1, p2;
    if (element.isVertical) {
        p1 = ImVec2(pos.x, pos.y - ELEMENT_HALF_LEN);
        p2 = ImVec2(pos.x, pos.y + ELEMENT_HALF_LEN);
    } else {
        p1 = ImVec2(pos.x - ELEMENT_HALF_LEN, pos.y);
        p2 = ImVec2(pos.x + ELEMENT_HALF_LEN, pos.y);
    }

    // Draw pin stubs and connection squares
    if (element.type != GROUND) {
        ImVec2 p_start, p_end;
        if (element.isVertical) {
            p_start = ImVec2(p1.x, p1.y - PIN_STUB);
            p_end = ImVec2(p2.x, p2.y + PIN_STUB);
        } else {
            p_start = ImVec2(p1.x - PIN_STUB, p1.y);
            p_end = ImVec2(p2.x + PIN_STUB, p2.y);
        }
        drawList->AddLine(p_start, p1, color, thickness);
        drawList->AddLine(p2, p_end, color, thickness);
        drawList->AddRectFilled(ImVec2(p_start.x - 2, p_start.y - 2), ImVec2(p_start.x + 2, p_start.y + 2), color);
        drawList->AddRectFilled(ImVec2(p_end.x - 2, p_end.y - 2), ImVec2(p_end.x + 2, p_end.y + 2), color);
    }

    // Draw the component's body based on its type and orientation
    switch (element.type) {
        case RESISTOR: {
            if (element.isVertical) {
                ImVec2 points[] = { p1, ImVec2(p1.x - 5, p1.y + 10), ImVec2(p1.x + 5, p1.y + 20),
                                    ImVec2(p1.x - 5, p1.y + 30), ImVec2(p1.x + 5, p1.y + 40),
                                    ImVec2(p1.x - 5, p1.y + 50), p2 };
                drawList->AddPolyline(points, 7, color, ImDrawFlags_None, thickness);
            } else {
                ImVec2 points[] = { p1, ImVec2(p1.x + 10, p1.y - 5), ImVec2(p1.x + 20, p1.y + 5),
                                    ImVec2(p1.x + 30, p1.y - 5), ImVec2(p1.x + 40, p1.y + 5),
                                    ImVec2(p1.x + 50, p1.y - 5), p2 };
                drawList->AddPolyline(points, 7, color, ImDrawFlags_None, thickness);
            }
            break;
        }

        case CAPACITOR: {
            if (element.isVertical) {
                drawList->AddLine(p1, ImVec2(pos.x, pos.y - 3), color, thickness);
                drawList->AddLine(ImVec2(pos.x - 10, pos.y - 3), ImVec2(pos.x + 10, pos.y - 3), color, thickness);
                drawList->AddLine(ImVec2(pos.x - 10, pos.y + 3), ImVec2(pos.x + 10, pos.y + 3), color, thickness);
                drawList->AddLine(ImVec2(pos.x, pos.y + 3), p2, color, thickness);
            } else {
                drawList->AddLine(p1, ImVec2(pos.x - 3, pos.y), color, thickness);
                drawList->AddLine(ImVec2(pos.x - 3, pos.y - 10), ImVec2(pos.x - 3, pos.y + 10), color, thickness);
                drawList->AddLine(ImVec2(pos.x + 3, pos.y - 10), ImVec2(pos.x + 3, pos.y + 10), color, thickness);
                drawList->AddLine(ImVec2(pos.x + 3, pos.y), p2, color, thickness);
            }
            break;
        }

        case INDUCTOR: {
            const int segments = 12;
            const float radius = 5.0f;
            if (element.isVertical) {
                drawList->AddLine(p1, ImVec2(p1.x, p1.y + 15), color, thickness);
                drawList->PathClear();
                for (int i = 0; i < 4; ++i) {
                    drawList->PathArcTo(ImVec2(p1.x, p1.y + 15 + (i * 2 + 1) * radius), radius, -M_PI / 2.0f, M_PI / 2.0f, segments);
                }
                drawList->PathStroke(color, ImDrawFlags_None, thickness);
                drawList->AddLine(ImVec2(p2.x, p2.y - 15), p2, color, thickness);
            } else {
                drawList->AddLine(p1, ImVec2(p1.x + 15, p1.y), color, thickness);
                drawList->PathClear();
                for (int i = 0; i < 4; ++i) {
                    drawList->PathArcTo(ImVec2(p1.x + 15 + (i * 2 + 1) * radius, p1.y), radius, M_PI, 2.0f * M_PI, segments);
                }
                drawList->PathStroke(color, ImDrawFlags_None, thickness);
                drawList->AddLine(ImVec2(p2.x - 15, p2.y), p2, color, thickness);
            }
            break;
        }

        case VSOURCE: {
            drawList->AddCircle(pos, 15.0f, color, 0, thickness);
            if (element.isVertical) { // Top is positive (+)
                drawList->AddLine(ImVec2(pos.x - 3, pos.y - 7), ImVec2(pos.x + 3, pos.y - 7), color, thickness);
                drawList->AddLine(ImVec2(pos.x, pos.y - 10), ImVec2(pos.x, pos.y - 4), color, thickness);
                drawList->AddLine(ImVec2(pos.x - 3, pos.y + 7), ImVec2(pos.x + 3, pos.y + 7), color, thickness);
            } else { // Left is positive (+)
                drawList->AddLine(ImVec2(pos.x - 10, pos.y), ImVec2(pos.x - 4, pos.y), color, thickness);
                drawList->AddLine(ImVec2(pos.x - 7, pos.y - 3), ImVec2(pos.x - 7, pos.y + 3), color, thickness);
                drawList->AddLine(ImVec2(pos.x + 4, pos.y), ImVec2(pos.x + 10, pos.y), color, thickness);
            }
            break;
        }

        case CSOURCE: {
            drawList->AddCircle(pos, 15.0f, color, 0, thickness);
            if (element.isVertical) { // Arrow points from top (p1) to bottom (p2)
                drawList->AddLine(ImVec2(pos.x, pos.y - 10), ImVec2(pos.x, pos.y + 10), color, thickness);
                drawList->AddLine(ImVec2(pos.x, pos.y + 10), ImVec2(pos.x - 4, pos.y + 5), color, thickness);
                drawList->AddLine(ImVec2(pos.x, pos.y + 10), ImVec2(pos.x + 4, pos.y + 5), color, thickness);
            } else { // Arrow points from left (p1) to right (p2)
                drawList->AddLine(ImVec2(pos.x - 10, pos.y), ImVec2(pos.x + 10, pos.y), color, thickness);
                drawList->AddLine(ImVec2(pos.x + 10, pos.y), ImVec2(pos.x + 5, pos.y - 4), color, thickness);
                drawList->AddLine(ImVec2(pos.x + 10, pos.y), ImVec2(pos.x + 5, pos.y + 4), color, thickness);
            }
            break;
        }

        case DIODE: {
            if (element.isVertical) { // Points down
                ImVec2 tri_p1 = ImVec2(pos.x - 8, p1.y + 15);
                ImVec2 tri_p2 = ImVec2(pos.x + 8, p1.y + 15);
                ImVec2 tri_p3 = ImVec2(pos.x, p2.y - 15);
                drawList->AddTriangle(tri_p1, tri_p2, tri_p3, color, thickness);
                drawList->AddLine(ImVec2(pos.x - 8, p2.y - 15), ImVec2(pos.x + 8, p2.y - 15), color, thickness);
            } else { // Points right
                ImVec2 tri_p1 = ImVec2(p1.x + 15, pos.y - 8);
                ImVec2 tri_p2 = ImVec2(p1.x + 15, pos.y + 8);
                ImVec2 tri_p3 = ImVec2(p2.x - 15, pos.y);
                drawList->AddTriangle(tri_p1, tri_p2, tri_p3, color, thickness);
                drawList->AddLine(ImVec2(p2.x - 15, pos.y - 8), ImVec2(p2.x - 15, pos.y + 8), color, thickness);
            }
            break;
        }

        case VSIN: {
            drawList->AddCircle(pos, 15.0f, color, 0, thickness);
            // Draw a small sine wave inside
            drawList->PathLineTo(ImVec2(pos.x - 10, pos.y));
            drawList->PathBezierCubicCurveTo(ImVec2(pos.x - 5, pos.y - 10), ImVec2(pos.x + 5, pos.y + 10), ImVec2(pos.x + 10, pos.y), 10);
            drawList->PathStroke(color, ImDrawFlags_None, thickness);
            break;
        }

        case CSIN: {
            drawList->AddCircle(pos, 15.0f, color, 0, thickness);
//            if (element.isVertical) { /* ... */ } else { /* ... */ }
//            // Add a small sine wave to distinguish it
            drawList->PathLineTo(ImVec2(pos.x - 7, pos.y - 7));
            drawList->PathBezierCubicCurveTo(ImVec2(pos.x - 3, pos.y - 11), ImVec2(pos.x + 3, pos.y - 3), ImVec2(pos.x + 7, pos.y - 7), 10);
            drawList->PathStroke(color, ImDrawFlags_None, thickness);
            break;
        }

        case GROUND: {
            drawList->AddLine(pos, ImVec2(pos.x, pos.y + 10), color, thickness);
            drawList->AddLine(ImVec2(pos.x - 10, pos.y + 10), ImVec2(pos.x + 10, pos.y + 10), color, thickness);
            drawList->AddLine(ImVec2(pos.x - 6, pos.y + 14), ImVec2(pos.x + 6, pos.y + 14), color, thickness);
            drawList->AddLine(ImVec2(pos.x - 3, pos.y + 18), ImVec2(pos.x + 3, pos.y + 18), color, thickness);
            break;
        }

        default:
            drawList->AddLine(p1, p2, color, thickness);
            break;
    }


    if (!element.name.empty() || element.value != 0 || element.type == VSIN) {
        char display_str[128];

        if (element.type == VSIN) {
            // Format each part individually using our helper function
            std::string offset_str = formatValueWithSI(element.value);
            std::string amp_str = formatValueWithSI(element.amp);
            std::string freq_str = formatValueWithSI(element.freq);

            // Combine the pre-formatted strings
            snprintf(display_str, 128, "SINE(%s %s %s)", offset_str.c_str(), amp_str.c_str(), freq_str.c_str());
        } else if (element.type == CSIN) {
            // Special formatting for AC Current Sources
            std::string offset_str = formatValueWithSI(element.value);
            std::string amp_str = formatValueWithSI(element.amp);
            std::string freq_str = formatValueWithSI(element.freq);
            snprintf(display_str, 128, "SINE(%s %s %s)", offset_str.c_str(), amp_str.c_str(), freq_str.c_str());
        } else {
            // Regular formatting for all other components
            snprintf(display_str, 128, "%s", formatValueWithSI(element.value).c_str());
        }

        ImVec2 name_offset = element.isVertical ? ImVec2(10, -15) : ImVec2(-15, -30);
        ImVec2 value_offset = element.isVertical ? ImVec2(10, 5) : ImVec2(-15, 15);

        drawList->AddText(ImVec2(pos.x + name_offset.x, pos.y + name_offset.y), color, element.name.c_str());
        drawList->AddText(ImVec2(pos.x + value_offset.x, pos.y + value_offset.y), color, display_str);
    }
}


void DrawProbe(ImDrawList* drawList, const PlacedProbe& probe, const ImVec2& canvas_p0) {
    const ImU32 color = IM_COL32(23, 107, 135, 255); // A nice blue color
    const ImU32 text_color = IM_COL32(0, 0, 0, 255);
    const float thickness = 1.5f;

    ImVec2 pos = ImVec2(canvas_p0.x + probe.pos.x, canvas_p0.y + probe.pos.y);
    drawList->AddCircle(pos, 15.0f, color, 0, thickness);
    drawList->AddText(ImVec2(pos.x - 4, pos.y - 8), color, "V");

    // Format the measured voltage text
    char buffer[32];
    snprintf(buffer, 32, "%.3f V", probe.measured_voltage);
    drawList->AddText(ImVec2(pos.x + 20, pos.y - 8), text_color, buffer);
}


// Snaps a position to the nearest grid point
ImVec2 SnapToGrid(const ImVec2& pos, float grid_step) {
    float snapped_x = roundf(pos.x / grid_step) * grid_step;
    float snapped_y = roundf(pos.y / grid_step) * grid_step;
    return ImVec2(snapped_x, snapped_y);
}


void UpdateAndDrawDebugInfo(ImDrawList* drawList, const ImVec2& canvas_p0) {
    if (!g_showDebugInfo) return;

    // This function re-runs the connection analysis to get up-to-date info.
    // NOTE: This is for debugging and is not performance-optimized.

    // 1. Reset and find all connectors
    std::vector<Connector> temp_connectors;
    std::map<int, std::vector<int>> temp_graph;
    int temp_counter = 0;

    auto add_conn_debug = [&](ImVec2 pos) -> int {
        for (const auto& c : temp_connectors) {
            float dx = c.pos.x - pos.x; float dy = c.pos.y - pos.y;
            if (dx * dx + dy * dy < 4.0f) return c.id;
        }
        temp_connectors.push_back({temp_counter, pos});
        return temp_counter++;
    };

    for (const auto& e : g_placedElements) {
        if (e.type == GROUND) {
            add_conn_debug(e.pos);
        } else {
            if (e.isVertical) {
                add_conn_debug(ImVec2(e.pos.x, e.pos.y - PIN_OFFSET));
                add_conn_debug(ImVec2(e.pos.x, e.pos.y + PIN_OFFSET));
            } else {
                add_conn_debug(ImVec2(e.pos.x - PIN_OFFSET, e.pos.y));
                add_conn_debug(ImVec2(e.pos.x + PIN_OFFSET, e.pos.y));
            }
        }
    }

    for (const auto& w : g_wires) {
        int id1 = add_conn_debug(w.p1); int id2 = add_conn_debug(w.p2);
        temp_graph[id1].push_back(id2); temp_graph[id2].push_back(id1);
    }

    // 2. Run DSU to find node groups
    std::map<int, int> parent;
    std::function<int(int)> find = [&](int i) -> int {
        if (parent.find(i) == parent.end() || parent[i] == i) return i;
        return parent[i] = find(parent[i]);
    };
    auto unite = [&](int i, int j) {
        int root_i = find(i); int root_j = find(j);
        if (root_i != root_j) parent[root_i] = root_j;
    };
    for (const auto& c : temp_connectors) parent[c.id] = c.id;
    for (const auto& pair : temp_graph) {
        for (int neighbor : pair.second) unite(pair.first, neighbor);
    }

    // 3. Find ground root ID
    int ground_root_id = -1;
    for (const auto& e : g_placedElements) {
        if (e.type == GROUND) { ground_root_id = find(add_conn_debug(e.pos)); break; }
    }

    // 4. Assign names and draw the debug info
    std::map<int, std::string> root_to_name;
    int node_idx_counter = 1;
    for (const auto& c : temp_connectors) {
        int root = find(c.id);
        if (root_to_name.find(root) == root_to_name.end()) {
            if (root == ground_root_id) root_to_name[root] = "0 (GND)";
            else root_to_name[root] = "n" + std::to_string(node_idx_counter++);
        }
        std::string node_name = root_to_name[root];
        ImVec2 pos_abs = ImVec2(canvas_p0.x + c.pos.x, canvas_p0.y + c.pos.y);
        drawList->AddCircleFilled(pos_abs, 4.0f, IM_COL32(255, 0, 0, 200));
        drawList->AddText(ImVec2(pos_abs.x + 5, pos_abs.y - 5), IM_COL32(255, 0, 0, 255), node_name.c_str());
    }
}


void RenderCanvas() {
    ImGui::BeginChild("Canvas", ImVec2(0, 0), true, ImGuiWindowFlags_NoMove);

    ImDrawList* drawList = ImGui::GetWindowDrawList();
    const ImVec2 canvas_p0 = ImGui::GetCursorScreenPos();
    const ImVec2 canvas_sz = ImGui::GetContentRegionAvail();
    const float GRID_STEP = 10.0f;

    // Draw grid
    for (float x = 0; x < canvas_sz.x; x += GRID_STEP) {
        for (float y = 0; y < canvas_sz.y; y += GRID_STEP) {
            drawList->AddCircleFilled(ImVec2(canvas_p0.x + x, canvas_p0.y + y), 1.0f, IM_COL32(200, 200, 200, 255));
        }
    }

    if (ImGui::IsWindowHovered()) {
        // Handle Rotation Key Press for previews
        if (g_currentTool != NONE && g_currentTool != WIRE && g_currentTool != VOLTMETER && ImGui::IsKeyPressed(ImGuiKey_R)) {
            g_previewIsVertical = !g_previewIsVertical;
        }

        const ImVec2 mouse_pos_in_canvas = ImVec2(ImGui::GetMousePos().x - canvas_p0.x, ImGui::GetMousePos().y - canvas_p0.y);
        const ImVec2 snapped_pos = SnapToGrid(mouse_pos_in_canvas, GRID_STEP);

        // --- Placement Logic ---
        if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
            if (g_currentTool == VOLTMETER) {
                buildCircuit();
                int closest_connector_id = -1;
                float min_dist_sq = 100.0f;
                for (const auto& conn : g_allConnectors) {
                    float dx = conn.pos.x - snapped_pos.x;
                    float dy = conn.pos.y - snapped_pos.y;
                    float dist_sq = dx * dx + dy * dy;
                    if (dist_sq < min_dist_sq) {
                        min_dist_sq = dist_sq;
                        closest_connector_id = conn.id;
                    }
                }
                if (closest_connector_id != -1) {
                    try {
                        g_circuit.computeDCOP(g_dcVoltages, g_dcCurrents);
                        std::string node_name = g_connectorToNodeName.at(closest_connector_id);
                        double voltage = (node_name == g_circuit.groundName) ? 0.0 : g_dcVoltages.at(node_name);
                        g_probes.push_back({snapped_pos, node_name, voltage});
                    } catch (const std::exception& e) {
                        g_simulationError = e.what();
                        g_showDCResults = true;
                    }
                }
                g_currentTool = NONE;
            } else if (g_currentTool == WIRE) {
                if (!g_isPlacingWire) {
                    g_isPlacingWire = true;
                    g_wireStartPos = snapped_pos;
                } else {
                    g_wires.push_back({g_wireStartPos, snapped_pos});
                    g_isPlacingWire = false;
                    g_currentTool = NONE;
                }
            } else if (g_currentTool != NONE) {
                std::string name;
                double value = 0, amp = 0, freq = 0;

                // Check if the tool is our new AC Source
                if (g_currentTool == VSIN) {
                    name = g_vacNameBuffer;
                    value = parseNumber(g_vacDcOffsetText); // DC offset goes into the main 'value' field
                    amp = parseNumber(g_vacAmpText);
                    freq = parseNumber(g_vacFreqText);
                } else if (g_currentTool == CSIN) { // <<< ADD THIS BLOCK
                    name = g_iacNameBuffer;
                    value = parseNumber(g_iacDcOffsetText);
                    amp = parseNumber(g_iacAmpText);
                    freq = parseNumber(g_iacFreqText);
                } else {
                    // Otherwise, use the regular component dialog buffers
                    name = g_componentNameBuffer;
                    value = parseNumber(g_componentValueBuffer);
                }

                if (g_currentTool == GROUND) {
                    name = "GND";
                    value = 0.0;
                }

                // Add the new element to our list with all its properties
                g_placedElements.push_back({g_currentTool, snapped_pos, name, value, g_previewIsVertical, amp, freq});

                g_currentTool = NONE;
                g_previewIsVertical = false;
            }
        }

        // --- Preview Logic ---
        if (g_currentTool != NONE) {
            if (g_currentTool == WIRE) {
                if (g_isPlacingWire) {
                    ImVec2 start_abs = ImVec2(canvas_p0.x + g_wireStartPos.x, canvas_p0.y + g_wireStartPos.y);
                    ImVec2 end_abs = ImVec2(canvas_p0.x + snapped_pos.x, canvas_p0.y + snapped_pos.y);
                    drawList->AddLine(start_abs, end_abs, IM_COL32(0, 0, 255, 255), 1.5f);
                }
            } else if (g_currentTool == VOLTMETER) {
                drawList->AddCircle(ImVec2(canvas_p0.x + snapped_pos.x, canvas_p0.y + snapped_pos.y), 15.0f, IM_COL32(23, 107, 135, 255), 0, 1.5f);
                drawList->AddText(ImVec2(canvas_p0.x + snapped_pos.x - 4, canvas_p0.y + snapped_pos.y - 8), IM_COL32(23, 107, 135, 255), "V");
            } else { // Component preview now includes VSIN
                PlacedElement preview_element = { g_currentTool, snapped_pos, "", 0.0, g_previewIsVertical };
                DrawComponent(drawList, preview_element, canvas_p0);
            }
        }
    }

    // --- Drawing Placed Items ---
    for (const auto& element : g_placedElements) DrawComponent(drawList, element, canvas_p0);
    for (const auto& wire : g_wires) {
        ImVec2 p1 = ImVec2(canvas_p0.x + wire.p1.x, canvas_p0.y + wire.p1.y);
        ImVec2 p2 = ImVec2(canvas_p0.x + wire.p2.x, canvas_p0.y + wire.p2.y);
        drawList->AddLine(p1, p2, IM_COL32(0, 0, 0, 255), 1.5f);
        drawList->AddRectFilled(ImVec2(p1.x - 2, p1.y - 2), ImVec2(p1.x + 2, p1.y + 2), IM_COL32(0, 0, 0, 255));
        drawList->AddRectFilled(ImVec2(p2.x - 2, p2.y - 2), ImVec2(p2.x + 2, p2.y + 2), IM_COL32(0, 0, 0, 255));
    }
    for (const auto& probe : g_probes) DrawProbe(drawList, probe, canvas_p0);

    UpdateAndDrawDebugInfo(drawList, canvas_p0);

    ImGui::EndChild();
}


// TODO: MAIN LOOP
int main() {
    // 1. Initialize GLFW and create a window
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return -1;
    }

    // Request OpenGL 3.3 context
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow* window = glfwCreateWindow(1280, 720, "Opulator Spice (ImGui Version)", NULL, NULL);
    if (!window) {
        std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1); // Enable vsync

    // 2. Initialize GLEW
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW" << std::endl;
        return -1;
    }

    // 3. Setup Dear ImGui and ImPlot contexts
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard; // Enable Keyboard Controls

    // Setup Dear ImGui style
    ImGui::StyleColorsLight();

    // Setup Platform/Renderer backends
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");


    // 4. Main application loop
    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // --- Create a single, fullscreen, non-interactable window as our main surface ---
        const ImGuiViewport* main_viewport = ImGui::GetMainViewport();
        ImGui::SetNextWindowPos(main_viewport->Pos);
        ImGui::SetNextWindowSize(main_viewport->Size);

        ImGuiWindowFlags window_flags =
                ImGuiWindowFlags_MenuBar |
                ImGuiWindowFlags_NoDecoration |
                ImGuiWindowFlags_NoMove |
                ImGuiWindowFlags_NoBringToFrontOnFocus |
                ImGuiWindowFlags_NoNavFocus;

        ImGui::Begin("MainApp", nullptr, window_flags);

        // Call our UI rendering functions inside this main window
        RenderMenuBar();
        RenderToolbar();
        RenderComponentPopup();
        RenderVACPopup();
        RenderIACPopup();
        RenderRunPopup();
        RenderDCResultsWindow();
        RenderTransientPlotWindow();
        RenderACPlotWindow();
        RenderCanvas();

        ImGui::End();
        // --- End of main surface window ---


        // Rendering
        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
    }

    // 5. Cleanup
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImPlot::DestroyContext();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}