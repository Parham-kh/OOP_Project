// ===========================================
//========== ماییم و نوای بی نوایی ==========
// ========== بسم الله اگرحریف مایی ==========
// ===========================================

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
#include "save_load.h"
#include "recent_files.h"

// Add with your other includes
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

// ================ Global Application State ================
// Enum to represent the currently selected tool
enum ToolType { NONE, RESISTOR, CAPACITOR, INDUCTOR, VSOURCE, CSOURCE, GROUND, WIRE, VSIN, CSIN, VOLTMETER, DIODE, VPULSE, IPULSE, VDELTA, IDELTA, VCVS, VCCS, VPHASE, SUBCIRCUIT_TOOL, SUBCIRCUIT_INSTANCE };
ToolType g_currentTool = NONE;


struct PlacedProbe {
    ImVec2 pos;
    std::string node_name;
    double measured_voltage;
};
std::vector<PlacedProbe> g_probes;


std::vector<SubcircuitDefinition> g_subcircuitLibrary;


// Struct to hold information about a placed graphical element
struct PlacedElement {
    ToolType type;
    ImVec2 pos;
    std::string name;
    double value;
    bool isVertical = false;

    // SINE parameters
    double amp = 0.0;
    double freq = 0.0;

    // PULSE parameters
    double v_initial = 0.0;
    double v_on = 0.0;

    double t_delay = 0.0;
    double t_rise = 0.0;
    double t_fall = 0.0;
    double t_on = 0.0;
    double t_period = 0.0;

    double i_initial = 0.0;
    double i_on = 0.0;

    // DELTA parameters
    double t_pulse = 0.0;
    double area = 1.0;

    // Dependent Source parameters
    double gain = 1.0;
    char ctrlNodeP_name[64] = "";
    char ctrlNodeN_name[64] = "";
    int subcircuit_def_idx = -1;
    ImVec2 size = {0,0};

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
enum RunChoice { RUN_DC, RUN_AC, RUN_VT, RUN_IT, RUN_PHASE };
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

bool g_showVPulsePopup = false;
char g_vpulseNameBuffer[64] = "V3";
char g_vpulseVInitial[64] = "0";
char g_vpulseVOn[64] = "5";
char g_vpulseTDelay[64] = "1s";
char g_vpulseTRise[64] = "10n";
char g_vpulseTFall[64] = "10n";
char g_vpulseTOn[64] = "1s";
char g_vpulseTPeriod[64] = "2s";

bool g_showIPulsePopup = false;
char g_ipulseNameBuffer[64] = "I3";
char g_ipulseIInitial[64] = "0";
char g_ipulseIOn[64] = "1";
char g_ipulseTDelay[64] = "1s";
char g_ipulseTRise[64] = "10n";
char g_ipulseTFall[64] = "10n";
char g_ipulseTOn[64] = "1s";
char g_ipulseTPeriod[64] = "2s";

bool g_showVDeltaPopup = false;
char g_vdeltaNameBuffer[64] = "Vd1";
char g_vdeltaTPulse[64] = "1s";
char g_vdeltaArea[64] = "1";

bool g_showIDeltaPopup = false;
char g_ideltaNameBuffer[64] = "Id1";
char g_ideltaTPulse[64] = "1s";
char g_ideltaArea[64] = "1";

bool g_showDepSourcePopup = false;
char g_depNameBuffer[64] = "G1";
char g_depGainBuffer[64] = "1";
char g_depCtrlPBuffer[64] = "n1";
char g_depCtrlNBuffer[64] = "n2";

std::string g_depPopupError;

int g_selectedIndex = -1; // -1 means no element is selected
bool g_isMovingElement = false;
ImVec2 g_originalPos;


// --- State for Plot Math Channel ---
// A list of all available signal names for dropdowns
std::vector<const char*> g_plot_signal_names;
// Maps dropdown index back to the signal name and type (V or I)
std::map<int, std::pair<std::string, char>> g_plot_signal_map;
// Indices of the two signals selected for the math operation
int g_mathSignal1_idx = -1;
int g_mathSignal2_idx = -1;
// Booleans to control visibility of the calculated traces
bool g_showAddTrace = false;
bool g_showSubtractTrace = false;
// Storage for the calculated data
std::vector<double> g_addTrace;
std::vector<double> g_subtractTrace;


// --- State for Measurement Cursors ---
bool g_cursor1_active = false;
bool g_cursor2_active = false;
ImPlotPoint g_cursor1_pos; // From implot.h, has .x and .y
ImPlotPoint g_cursor2_pos;
int g_cursor1_series_idx = -1; // Which signal the cursor is attached to
int g_cursor2_series_idx = -1;

std::string g_currentFilePath;

int g_subcircuitSelectionStep = 0; // 0=none, 1=first port selected
std::string g_subcircuitNode1_name, g_subcircuitNode2_name;
ImVec2 g_subcircuitNode1_pos, g_subcircuitNode2_pos;
char g_subcircuitNameBuffer[64] = "MySub";
int g_placingSubcircuit_idx = -1;

bool g_showVPhasePopup = false;
char g_vphaseNameBuffer[64] = "Vp1";
char g_vphaseAmpBuffer[64] = "1";
char g_vphaseFreqBuffer[64] = "1k";
char g_vphasePhaseBuffer[64] = "0";

bool g_showPhasePlot = false;
char g_phaseBaseFreqText[64] = "1k";
char g_phaseStartText[64] = "0";
char g_phaseStopText[64] = "360"; // Degrees are more intuitive for UI
char g_phaseNumStepsText[64] = "100";
char g_phaseProbeText[64] = "n1";
std::vector<double> g_phase_angles; // In degrees
std::vector<double> g_phase_mags;
std::string g_phase_probe_node;

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
    enum class Type { DC, SINE, PULSE, DELTA, VCVS };
    Type sourceType = Type::DC;

    // SINE parameters
    double dcOffset = 0.0;
    double amplitude = 0.0;
    double frequency = 0.0;
    double phase = 0.0;

    // PULSE parameters
    double v_initial = 0.0;
    double v_on = 0.0;
    double t_delay = 0.0;
    double t_rise = 0.0;
    double t_fall = 0.0;
    double t_on = 0.0;
    double t_period = 0.0;

    // DELTA parameters
    double t_pulse = 0.0;
    double area = 1.0;

    // Dependent source parameters
    double gain = 1.0;
    Node* ctrlNodeP = nullptr;
    Node* ctrlNodeN = nullptr;

    // Constructors
    VoltageSource(const std::string &elemName, Node* n1, Node* n2, double val)
            : Element(elemName, n1, n2, val) {
        sourceType = Type::DC;
        dcOffset = val;
    }

    VoltageSource(const std::string &elemName, Node* n1, Node* n2, double offset, double amp, double freq, double ph = 0.0)
            : Element(elemName, n1, n2, offset), dcOffset(offset), amplitude(amp), frequency(freq), phase(ph) {
        sourceType = Type::SINE;
    }

    std::string getType() const override { return "VoltageSource"; }

    // --- THIS IS THE CORRECTED FUNCTION SIGNATURE ---
    void setValueAtTime(double t, double step) {
        switch (sourceType) {
            case Type::SINE:
                this->value = dcOffset + amplitude * sin(2 * M_PI * frequency * t);
                break;
            case Type::PULSE: {
                if (t < t_delay) { this->value = v_initial; return; }
                double localTime = t - t_delay;
                double cycleTime = fmod(localTime, t_period);
                if (t_rise > 0 && cycleTime < t_rise) {
                    this->value = v_initial + (v_on - v_initial) * (cycleTime / t_rise);
                } else if (cycleTime < t_rise + t_on) {
                    this->value = v_on;
                } else if (t_fall > 0 && cycleTime < t_rise + t_on + t_fall) {
                    this->value = v_on - (v_on - v_initial) * ((cycleTime - t_rise - t_on) / t_fall);
                } else {
                    this->value = v_initial;
                }
                break;
            }
            case Type::DELTA:
                // This logic now works because it receives the 'step' value
                if (std::abs(t - t_pulse) < (step / 2.0)) {
                    this->value = area / step;
                } else {
                    this->value = 0.0;
                }
                break;
            case Type::DC:
            default:
                this->value = dcOffset;
                break;
        }
    }
};


class CurrentSource : public Element {
public:
    enum class Type { DC, SINE, PULSE, DELTA, VCCS };
    Type sourceType = Type::DC;

    // SINE/DC parameters
    double dcOffset = 0.0;
    double amplitude = 0.0;
    double frequency = 0.0;

    // PULSE parameters
    double i_initial = 0.0;
    double i_on = 0.0;
    double t_delay = 0.0;
    double t_rise = 0.0;
    double t_fall = 0.0;
    double t_on = 0.0;
    double t_period = 0.0;

    // DELTA parameters
    double t_pulse = 0.0;
    double area = 1.0;

    // Dependent source parameters
    double gain = 1.0;
    Node* ctrlNodeP = nullptr;
    Node* ctrlNodeN = nullptr;


    CurrentSource(const std::string &elemName, Node* a, Node* b, double offset, double amp = 0, double freq = 0)
            : Element(elemName, a, b, offset), dcOffset(offset), amplitude(amp), frequency(freq) {
        if (amp != 0 && freq != 0) sourceType = Type::SINE;
        else sourceType = Type::DC;
    }

    std::string getType() const override { return "CurrentSource"; }

    void setValueAtTime(double t, double step) {
        switch (sourceType) {
            case Type::SINE:
                this->value = dcOffset + amplitude * sin(2 * M_PI * frequency * t);
                break;
            case Type::PULSE: {
                if (t < t_delay) { this->value = i_initial; return; }
                double localTime = t - t_delay;
                double cycleTime = fmod(localTime, t_period);
                if (t_rise > 0 && cycleTime < t_rise) {
                    this->value = i_initial + (i_on - i_initial) * (cycleTime / t_rise);
                } else if (cycleTime < t_rise + t_on) {
                    this->value = i_on;
                } else if (t_fall > 0 && cycleTime < t_rise + t_on + t_fall) {
                    this->value = i_on - (i_on - i_initial) * ((cycleTime - t_rise - t_on) / t_fall);
                } else {
                    this->value = i_initial;
                }
            } break;
            case Type::DELTA:
                if (std::abs(t - t_pulse) < (step / 2.0)) {
                    this->value = area / step;
                } else {
                    this->value = 0.0;
                }
                break;
            case Type::DC:
            default:
                this->value = dcOffset;
                break;
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
                    cs->setValueAtTime(t, step);
                    double I_inst = cs->getValue();
                    if (pi != -1) I[pi] -= I_inst;
                    if (ni != -1) I[ni] += I_inst;
                }
            }

            for (int k=0; k<m; ++k) {
                auto vs = voltageSources[k];
                Node* pn = vs->getNode1();
                Node* nn = vs->getNode2();
                int pi = (pn->getName() != this->groundName) ? nodeIndex.at(pn->getName()) : -1;
                int ni = (nn->getName() != this->groundName) ? nodeIndex.at(nn->getName()) : -1;

                vs->setValueAtTime(t, step);
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

                if (auto vccs = dynamic_cast<CurrentSource*>(elem); vccs && vccs->sourceType == CurrentSource::Type::VCCS) {
                    double g = vccs->gain;
                    int idx_p = (vccs->getNode1()->getName() != groundName) ? node_idx.at(vccs->getNode1()->getName()) : -1;
                    int idx_n = (vccs->getNode2()->getName() != groundName) ? node_idx.at(vccs->getNode2()->getName()) : -1;
                    int idx_cp = (vccs->ctrlNodeP->getName() != groundName) ? node_idx.at(vccs->ctrlNodeP->getName()) : -1;
                    int idx_cn = (vccs->ctrlNodeN->getName() != groundName) ? node_idx.at(vccs->ctrlNodeN->getName()) : -1;

                    if (idx_p != -1 && idx_cp != -1) A[idx_p][idx_cp] += g;
                    if (idx_p != -1 && idx_cn != -1) A[idx_p][idx_cn] -= g;
                    if (idx_n != -1 && idx_cp != -1) A[idx_n][idx_cp] -= g;
                    if (idx_n != -1 && idx_cn != -1) A[idx_n][idx_cn] += g;
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

            // --- D. Stamp Voltage Sources, ON Diodes, and VCVS using MNA ---
            for (int k = 0; k < m; ++k) {
                const auto* elem = v_sources_and_on_diodes[k];
                std::string name_a = elem->getNode1()->getName();
                std::string name_c = elem->getNode2()->getName();
                int current_idx = n + k;

                int idx_a = (name_a != groundName) ? node_idx.at(name_a) : -1;
                int idx_c = (name_c != groundName) ? node_idx.at(name_c) : -1;

                if (auto vcvs = dynamic_cast<const VoltageSource*>(elem); vcvs && vcvs->sourceType == VoltageSource::Type::VCVS) {
                    // Add current contribution to KCL equations
                    if (idx_a != -1) A[idx_a][current_idx] += 1.0; // <<< FIX
                    if (idx_c != -1) A[idx_c][current_idx] -= 1.0; // <<< FIX

                    // Add voltage constraint equation: Va - Vc - gain*Vcp + gain*Vcn = 0
                    int idx_cp = (vcvs->ctrlNodeP->getName() != groundName) ? node_idx.at(vcvs->ctrlNodeP->getName()) : -1;
                    int idx_cn = (vcvs->ctrlNodeN->getName() != groundName) ? node_idx.at(vcvs->ctrlNodeN->getName()) : -1;

                    if (idx_a != -1) A[current_idx][idx_a] += 1.0;
                    if (idx_c != -1) A[current_idx][idx_c] -= 1.0;
                    if (idx_cp != -1) A[current_idx][idx_cp] -= vcvs->gain;
                    if (idx_cn != -1) A[current_idx][idx_cn] += vcvs->gain;

                    b[current_idx] = 0.0;
                } else {
                    // --- This is the logic for standard V-sources and ON-Diodes ---
                    double V = elem->getValue();
                    if (idx_a != -1) A[idx_a][current_idx] += 1.0;
                    if (idx_c != -1) A[idx_c][current_idx] -= 1.0;

                    if (idx_a != -1) A[current_idx][idx_a] += 1.0;
                    if (idx_c != -1) A[current_idx][idx_c] -= 1.0;

                    b[current_idx] = V;
                }
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


// Returns the standard path for the library file
std::string subcircuit_library_path() {
#ifdef _WIN32
    const char* appdata = std::getenv("APPDATA");
    if (appdata && *appdata) {
        return std::string(appdata) + "\\Opulator\\library.json";
    }
#endif
    return "library.json"; // Fallback
}

// Saves the current in-memory library to the file
void SaveSubcircuitLibrary() {
    std::string path = subcircuit_library_path();
    ensure_directory_exists(path); // From recent_files.h

    SubcircuitLibraryStore store;
    store.library = g_subcircuitLibrary; // Copy from our global variable

    std::ofstream os(path);
    if (!os) return;
    cereal::JSONOutputArchive ar(os);
    ar(cereal::make_nvp("subcircuit_library", store));
}

// Loads the library from the file into memory
void LoadSubcircuitLibrary() {
    std::string path = subcircuit_library_path();
    std::ifstream is(path);
    if (!is) return; // File doesn't exist yet, that's fine

    try {
        cereal::JSONInputArchive ar(is);
        SubcircuitLibraryStore store;
        ar(cereal::make_nvp("subcircuit_library", store));
        g_subcircuitLibrary = store.library; // Load into our global variable
    } catch (const std::exception& e) {
        std::cerr << "Error loading subcircuit library: " << e.what() << std::endl;
    }
}


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
        b[current_idx] = std::polar(vs->amplitude, vs->phase);
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



void buildCircuit(bool require_ground = true) {
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

    if (require_ground && ground_root_id == -1) {
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

//    // 6. Build the logical circuit object
//    g_circuit.setGroundNode("0");
//    for (const auto& pair : root_to_name) {
//        g_circuit.getOrCreateNode(pair.second);
//    }
//
//    for (const auto& e : g_placedElements) {
//        if (e.type == GROUND || e.type == WIRE || e.type == VOLTMETER) continue;
//
//        int id1, id2;
//        if (e.isVertical) {
//            id1 = addConnector(ImVec2(e.pos.x, e.pos.y - PIN_OFFSET));
//            id2 = addConnector(ImVec2(e.pos.x, e.pos.y + PIN_OFFSET));
//        } else {
//            id1 = addConnector(ImVec2(e.pos.x - PIN_OFFSET, e.pos.y));
//            id2 = addConnector(ImVec2(e.pos.x + PIN_OFFSET, e.pos.y));
//        }
//
//        Node* n1 = g_circuit.getOrCreateNode(g_connectorToNodeName.at(id1));
//        Node* n2 = g_circuit.getOrCreateNode(g_connectorToNodeName.at(id2));
//
//        switch (e.type) {
//            case RESISTOR: g_circuit.addElement(new Resistor(e.name, n1, n2, e.value)); break;
//            case CAPACITOR: g_circuit.addElement(new Capacitor(e.name, n1, n2, e.value)); break;
//            case INDUCTOR: g_circuit.addElement(new Inductor(e.name, n1, n2, e.value)); break;
//            case VSOURCE: g_circuit.addElement(new VoltageSource(e.name, n1, n2, e.value)); break;
//            case CSOURCE: g_circuit.addElement(new CurrentSource(e.name, n1, n2, e.value)); break;
//            case DIODE: g_circuit.addElement(new Diode(e.name, n1, n2, e.value)); break;
//            case VSIN: g_circuit.addElement(new VoltageSource(e.name, n1, n2, e.value, e.amp, e.freq)); break;
//            case CSIN: g_circuit.addElement(new CurrentSource(e.name, n1, n2, e.value, e.amp, e.freq)); break;
//            case VPULSE: {
//                auto vs = new VoltageSource(e.name, n1, n2, 0.0); // Create a base VSource
//                vs->sourceType = VoltageSource::Type::PULSE;
//                // Copy all the pulse parameters from the PlacedElement to the simulation object
//                vs->v_initial = e.v_initial;
//                vs->v_on = e.v_on;
//                vs->t_delay = e.t_delay;
//                vs->t_rise = e.t_rise;
//                vs->t_fall = e.t_fall;
//                vs->t_on = e.t_on;
//                vs->t_period = e.t_period;
//                g_circuit.addElement(vs);
//                break;
//            }
//            case IPULSE: {
//                auto cs = new CurrentSource(e.name, n1, n2, e.i_initial);
//                cs->sourceType = CurrentSource::Type::PULSE;
//                cs->i_initial = e.i_initial;
//                cs->i_on      = e.i_on;
//                cs->t_delay   = e.t_delay;
//                cs->t_rise    = e.t_rise;
//                cs->t_fall    = e.t_fall;
//                cs->t_on      = e.t_on;
//                cs->t_period  = e.t_period;
//                g_circuit.addElement(cs);
//                break;
//            }
//            case VDELTA: {
//                auto vs = new VoltageSource(e.name, n1, n2, 0.0); // DC value is 0
//                vs->sourceType = VoltageSource::Type::DELTA;
//                vs->t_pulse = e.t_pulse;
//                vs->area = e.area;
//                g_circuit.addElement(vs);
//                break;
//            }
//            case VCVS: {
//                auto vs = new VoltageSource(e.name, n1, n2, 0.0);
//                vs->sourceType = VoltageSource::Type::VCVS;
//                vs->gain = e.gain;
//                vs->ctrlNodeP = g_circuit.getOrCreateNode(e.ctrlNodeP_name);
//                vs->ctrlNodeN = g_circuit.getOrCreateNode(e.ctrlNodeN_name);
//                g_circuit.addElement(vs);
//                break;
//            }
//            case VCCS: {
//                auto cs = new CurrentSource(e.name, n1, n2, 0.0);
//                cs->sourceType = CurrentSource::Type::VCCS;
//                cs->gain = e.gain;
//                cs->ctrlNodeP = g_circuit.getOrCreateNode(e.ctrlNodeP_name);
//                cs->ctrlNodeN = g_circuit.getOrCreateNode(e.ctrlNodeN_name);
//                g_circuit.addElement(cs);
//                break;
//            }
//            default: break;
//        }
//    }
//
//    // Add wires as ideal resistors
//    for (const auto& w : g_wires) {
//        int id1 = addConnector(w.p1);
//        int id2 = addConnector(w.p2);
//        std::string name1 = g_connectorToNodeName.at(id1);
//        std::string name2 = g_connectorToNodeName.at(id2);
//        if (name1 != name2) {
//            Node* n1 = g_circuit.getOrCreateNode(name1);
//            Node* n2 = g_circuit.getOrCreateNode(name2);
//            g_circuit.addElement(new Resistor("W", n1, n2, 1e-9));
//        }
//    }

    // 6. Build the logical circuit object
    g_circuit.setGroundNode("0");
    for (const auto& pair : root_to_name) {
        g_circuit.getOrCreateNode(pair.second);
    }

    // 7. Add all top-level elements and EXPAND subcircuits
    for (const auto& e : g_placedElements) {
        if (e.type == GROUND || e.type == WIRE || e.type == VOLTMETER) continue;

        // Find the main nodes the element is connected to on the canvas
        Node* n1 = nullptr;
        Node* n2 = nullptr;
        if (e.type != GROUND) {
            ImVec2 p1_pos = e.isVertical ? ImVec2(e.pos.x, e.pos.y - PIN_OFFSET) : ImVec2(e.pos.x - PIN_OFFSET, e.pos.y);
            ImVec2 p2_pos = e.isVertical ? ImVec2(e.pos.x, e.pos.y + PIN_OFFSET) : ImVec2(e.pos.x + PIN_OFFSET, e.pos.y);
            n1 = g_circuit.getOrCreateNode(g_connectorToNodeName.at(addConnector(p1_pos)));
            n2 = g_circuit.getOrCreateNode(g_connectorToNodeName.at(addConnector(p2_pos)));
        }

        if (e.type == SUBCIRCUIT_INSTANCE) {
            // --- THIS IS THE NEW EXPANSION LOGIC ---
            const SubcircuitDefinition& def = g_subcircuitLibrary.at(e.subcircuit_def_idx);
            std::map<int, Node*> internal_node_map; // Map internal subcircuit node IDs to main circuit Node*

            // A. Find the internal integer IDs for the named ports
            int port1_id = -1, port2_id = -1;
            for (const auto& snode : def.Snodes) {
                if (snode.name == def.port1_node_name) port1_id = snode.id;
                if (snode.name == def.port2_node_name) port2_id = snode.id;
            }

            // B. Map the external ports to the main circuit nodes
            if (port1_id != -1) internal_node_map[port1_id] = n1;
            if (port2_id != -1) internal_node_map[port2_id] = n2;

            // C. Instantiate the subcircuit's internal elements
            for (const auto& sub_elem : def.Selements) {
                Node *sub_n1 = nullptr, *sub_n2 = nullptr;

                // Remap the first pin's node
                if (internal_node_map.count(sub_elem.n1)) {
                    sub_n1 = internal_node_map.at(sub_elem.n1);
                } else {
                    std::string new_node_name = e.name + "." + def.Snodes.at(sub_elem.n1).name;
                    sub_n1 = g_circuit.getOrCreateNode(new_node_name);
                    internal_node_map[sub_elem.n1] = sub_n1;
                }
                // Remap the second pin's node
                if (internal_node_map.count(sub_elem.n2)) {
                    sub_n2 = internal_node_map.at(sub_elem.n2);
                } else {
                    std::string new_node_name = e.name + "." + def.Snodes.at(sub_elem.n2).name;
                    sub_n2 = g_circuit.getOrCreateNode(new_node_name);
                    internal_node_map[sub_elem.n2] = sub_n2;
                }

                // Create a new instance of the internal element with a unique name
                std::string new_elem_name = e.name + "." + sub_elem.name;

                // --- THIS IS THE NEW, COMPLETE SWITCH STATEMENT ---
                if (sub_elem.kind == "R") {
                    g_circuit.addElement(new Resistor(new_elem_name, sub_n1, sub_n2, sub_elem.value));
                } else if (sub_elem.kind == "C") {
                    g_circuit.addElement(new Capacitor(new_elem_name, sub_n1, sub_n2, sub_elem.value));
                } else if (sub_elem.kind == "L") {
                    g_circuit.addElement(new Inductor(new_elem_name, sub_n1, sub_n2, sub_elem.value));
                } else if (sub_elem.kind == "D") {
                    g_circuit.addElement(new Diode(new_elem_name, sub_n1, sub_n2, sub_elem.value));
                }
                    // Handle complex sources inside the subcircuit
                else if (sub_elem.kind == "V") {
                    auto vs = new VoltageSource(new_elem_name, sub_n1, sub_n2, sub_elem.value);
                    if (sub_elem.sourceType == "SINE") {
                        vs->sourceType = VoltageSource::Type::SINE;
                        vs->amplitude = sub_elem.amp;
                        vs->frequency = sub_elem.freq;
                    } // ... Add cases for PULSE, DELTA, VCVS ...
                    g_circuit.addElement(vs);
                } else if (sub_elem.kind == "I") {
                    auto cs = new CurrentSource(new_elem_name, sub_n1, sub_n2, sub_elem.value);
                    if (sub_elem.sourceType == "SINE") {
                        cs->sourceType = CurrentSource::Type::SINE;
                        cs->amplitude = sub_elem.amp;
                        cs->frequency = sub_elem.freq;
                    } // ... Add cases for PULSE, DELTA, VCCS ...
                    g_circuit.addElement(cs);
                }
            }

        } else {
            // --- This is the existing logic for all other primitive components ---
            switch (e.type) {
                case RESISTOR: g_circuit.addElement(new Resistor(e.name, n1, n2, e.value)); break;
                case CAPACITOR: g_circuit.addElement(new Capacitor(e.name, n1, n2, e.value)); break;
                case INDUCTOR: g_circuit.addElement(new Inductor(e.name, n1, n2, e.value)); break;
                case VSOURCE: g_circuit.addElement(new VoltageSource(e.name, n1, n2, e.value)); break;
                case CSOURCE: g_circuit.addElement(new CurrentSource(e.name, n1, n2, e.value)); break;
                case DIODE: g_circuit.addElement(new Diode(e.name, n1, n2, e.value)); break;
                case VSIN: g_circuit.addElement(new VoltageSource(e.name, n1, n2, e.value, e.amp, e.freq)); break;
                case CSIN: g_circuit.addElement(new CurrentSource(e.name, n1, n2, e.value, e.amp, e.freq)); break;
                case VPULSE: {
                    auto vs = new VoltageSource(e.name, n1, n2, 0.0); // Create a base VSource
                    vs->sourceType = VoltageSource::Type::PULSE;
                    // Copy all the pulse parameters from the PlacedElement to the simulation object
                    vs->v_initial = e.v_initial;
                    vs->v_on = e.v_on;
                    vs->t_delay = e.t_delay;
                    vs->t_rise = e.t_rise;
                    vs->t_fall = e.t_fall;
                    vs->t_on = e.t_on;
                    vs->t_period = e.t_period;
                    g_circuit.addElement(vs);
                    break;
                }
                case IPULSE: {
                    auto cs = new CurrentSource(e.name, n1, n2, e.i_initial);
                    cs->sourceType = CurrentSource::Type::PULSE;
                    cs->i_initial = e.i_initial;
                    cs->i_on      = e.i_on;
                    cs->t_delay   = e.t_delay;
                    cs->t_rise    = e.t_rise;
                    cs->t_fall    = e.t_fall;
                    cs->t_on      = e.t_on;
                    cs->t_period  = e.t_period;
                    g_circuit.addElement(cs);
                    break;
                }
                case VDELTA: {
                    auto vs = new VoltageSource(e.name, n1, n2, 0.0); // DC value is 0
                    vs->sourceType = VoltageSource::Type::DELTA;
                    vs->t_pulse = e.t_pulse;
                    vs->area = e.area;
                    g_circuit.addElement(vs);
                    break;
                }
                case IDELTA: {
                    auto cs = new CurrentSource(e.name, n1, n2, 0.0); // Create a CurrentSource, not a VoltageSource
                    cs->sourceType = CurrentSource::Type::DELTA; // Use the CurrentSource's enum
                    cs->t_pulse = e.t_pulse;
                    cs->area = e.area;
                    g_circuit.addElement(cs); // Add the correct 'cs' object
                    break;
                }
                case VPHASE: // VPHASE is a type of SINE source with phase
                    g_circuit.addElement(new VoltageSource(e.name, n1, n2, e.value, e.amp, e.freq, e.phase));
                    break;
                case VCVS: {
                    auto vs = new VoltageSource(e.name, n1, n2, 0.0);
                    vs->sourceType = VoltageSource::Type::VCVS;
                    vs->gain = e.gain;
                    vs->ctrlNodeP = g_circuit.getOrCreateNode(e.ctrlNodeP_name);
                    vs->ctrlNodeN = g_circuit.getOrCreateNode(e.ctrlNodeN_name);
                    g_circuit.addElement(vs);
                    break;
                }
                case VCCS: {
                    auto cs = new CurrentSource(e.name, n1, n2, 0.0);
                    cs->sourceType = CurrentSource::Type::VCCS;
                    cs->gain = e.gain;
                    cs->ctrlNodeP = g_circuit.getOrCreateNode(e.ctrlNodeP_name);
                    cs->ctrlNodeN = g_circuit.getOrCreateNode(e.ctrlNodeN_name);
                    g_circuit.addElement(cs);
                    break;
                }
                default: break;
            }
        }
    }

    // 8. Add top-level wires (unchanged)
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


// This function converts the visual schematic into a savable SCircuit object
SCircuit createSaveDataFromCanvas(bool require_ground = true) {
    SCircuit ckt_data;

    // 1. Perform connection analysis to find all logical nodes
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

    int ground_root_id = -1;
    for (const auto& e : g_placedElements) {
        if (e.type == GROUND) {
            int conn_id = addConnector(e.pos);
            ground_root_id = find(conn_id);
            break;
        }
    }

    if (require_ground && ground_root_id == -1) { // <<< CHANGED
        throw NoGroundException();
    }

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

    // 2. Create the list of SNodes for the save file
    std::map<std::string, int> node_name_to_id;
    int next_node_id = 0;
    for (const auto& pair : root_to_name) {
        std::string name = pair.second;
        node_name_to_id[name] = next_node_id;

        // Find a representative position for the node for saving
        ImVec2 node_pos = {0,0};
        for (const auto& c : g_allConnectors) {
            if (find(c.id) == pair.first) {
                node_pos = c.pos;
                break;
            }
        }
        ckt_data.Snodes.push_back({next_node_id, name, (int)node_pos.x, (int)node_pos.y});
        next_node_id++;
    }

    for (const auto& e : g_placedElements) {
        if (e.type == WIRE || e.type == VOLTMETER) continue;

        SElement se;
        se.name = e.name;
        se.x = (int)e.pos.x;
        se.y = (int)e.pos.y;
        se.isVertical = e.isVertical;
        se.value = e.value;

        // Find node IDs for the element's pins
        if (e.type != GROUND) {
            ImVec2 p1_pos = e.isVertical ? ImVec2(e.pos.x, e.pos.y - PIN_OFFSET) : ImVec2(e.pos.x - PIN_OFFSET, e.pos.y);
            ImVec2 p2_pos = e.isVertical ? ImVec2(e.pos.x, e.pos.y + PIN_OFFSET) : ImVec2(e.pos.x + PIN_OFFSET, e.pos.y);
            se.n1 = node_name_to_id.at(g_connectorToNodeName.at(addConnector(p1_pos)));
            se.n2 = node_name_to_id.at(g_connectorToNodeName.at(addConnector(p2_pos)));
        } else {
            se.n1 = node_name_to_id.at(g_connectorToNodeName.at(addConnector(e.pos)));
            se.n2 = se.n1;
        }

        // Map ToolType enum to a string for saving
        switch (e.type) {
            case RESISTOR:  se.kind = "R"; break;
            case CAPACITOR: se.kind = "C"; break;
            case INDUCTOR:  se.kind = "L"; break;
            case DIODE:     se.kind = "D"; break;
            case GROUND:    se.kind = "GND"; break;
            case VSOURCE: case VSIN: case VPULSE: case VDELTA: case VCVS:
                se.kind = "V"; break;
            case CSOURCE: case CSIN: case IPULSE: case IDELTA: case VCCS:
                se.kind = "I"; break;
        }

        // Save source-specific parameters
        if (e.type == VSOURCE || e.type == CSOURCE) { se.sourceType = "DC"; }
        else if (e.type == VSIN || e.type == CSIN) { se.sourceType = "SINE"; se.amp = e.amp; se.freq = e.freq; }
        else if (e.type == VPULSE) {
            se.sourceType = "PULSE";
            se.v_initial = e.v_initial; se.v_on = e.v_on;
            se.t_delay = e.t_delay; se.t_rise = e.t_rise; se.t_fall = e.t_fall; se.t_on = e.t_on; se.t_period = e.t_period;
        } else if (e.type == IPULSE) {
            se.sourceType = "PULSE";
            se.i_initial = e.i_initial; se.i_on = e.i_on;
            se.t_delay = e.t_delay; se.t_rise = e.t_rise; se.t_fall = e.t_fall; se.t_on = e.t_on; se.t_period = e.t_period;
        }
        else if (e.type == VDELTA || e.type == IDELTA) {
            se.sourceType = "DELTA";
            se.t_pulse = e.t_pulse; se.area = e.area;
        } else if (e.type == VCVS || e.type == VCCS) {
            se.sourceType = (e.type == VCVS) ? "VCVS" : "VCCS";
            se.gain = e.gain;
            se.ctrlNodeP_name = e.ctrlNodeP_name;
            se.ctrlNodeN_name = e.ctrlNodeN_name;
        }

        ckt_data.Selements.push_back(se);
    }

    for (const auto& w : g_wires) {
        SWire sw;
        sw.points.push_back({(int)w.p1.x, (int)w.p1.y});
        sw.points.push_back({(int)w.p2.x, (int)w.p2.y});
        ckt_data.Swires.push_back(sw);
    }

    return ckt_data;
}


void NewProject() {
    g_placedElements.clear();
    g_wires.clear();
    g_probes.clear();
    g_selectedIndex = -1;
    g_isMovingElement = false;
    g_currentFilePath = "";
    g_circuit.reset();
}


void loadCanvasFromSaveData(const SCircuit& ckt_data) {
    NewProject(); // Start with a clean slate

    // Helper map to convert saved strings back to ToolType enums
    std::map<std::string, ToolType> kind_to_tool;
    kind_to_tool["R"] = RESISTOR;
    kind_to_tool["C"] = CAPACITOR;
    kind_to_tool["L"] = INDUCTOR;
    kind_to_tool["D"] = DIODE;
    kind_to_tool["GND"] = GROUND;

    std::map<std::string, ToolType> source_kind_to_tool;
    source_kind_to_tool["DC"] = VSOURCE; // Default
    source_kind_to_tool["SINE"] = VSIN;
    source_kind_to_tool["PULSE"] = VPULSE;
    source_kind_to_tool["DELTA"] = VDELTA;
    source_kind_to_tool["VCVS"] = VCVS;
    source_kind_to_tool["VCCS"] = VCCS;


    for (const auto& se : ckt_data.Selements) {
        PlacedElement e{};
        e.name = se.name;
        e.pos = ImVec2((float)se.x, (float)se.y);
        e.isVertical = se.isVertical;
        e.value = se.value;

        if (se.kind == "V") {
            e.type = source_kind_to_tool.count(se.sourceType) ? source_kind_to_tool[se.sourceType] : VSOURCE;
            // Copy all source parameters from se to e
            e.amp = se.amp; e.freq = se.freq;
            e.v_initial = se.v_initial; e.v_on = se.v_on;
            e.t_delay = se.t_delay; e.t_rise = se.t_rise; e.t_fall = se.t_fall; e.t_on = se.t_on; e.t_period = se.t_period;
            e.t_pulse = se.t_pulse; e.area = se.area;
            e.gain = se.gain;
            strcpy_s(e.ctrlNodeP_name, se.ctrlNodeP_name.c_str());
            strcpy_s(e.ctrlNodeN_name, se.ctrlNodeN_name.c_str());
        }
        else if (se.kind == "I") {
            // A similar map would be needed for I-sources if they have different ToolTypes
            e.type = CSOURCE; // Simplified for now
        }
        else {
            e.type = kind_to_tool.at(se.kind);
        }

        g_placedElements.push_back(e);
    }

    for (const auto& sw : ckt_data.Swires) {
        if (sw.points.size() >= 2) {
            g_wires.push_back({ImVec2((float)sw.points[0].x, (float)sw.points[0].y),
                               ImVec2((float)sw.points[1].x, (float)sw.points[1].y)});
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


void RenderSubcircuitMenu();

// ================ UI Rendering Functions ================
void DoSaveAs() {
    auto path = pfd::save_file("Save Circuit As", "", { "JSON Files", "*.json" }).result();
    if (!path.empty()) {
        if (path.rfind(".json") == std::string::npos && path.rfind(".JSON") == std::string::npos) {
            path += ".json";
        }
        g_currentFilePath = path;
        SCircuit ckt_data = createSaveDataFromCanvas(false);
        saveProject(g_currentFilePath, ckt_data);
        recent_add(g_currentFilePath);
    }
}

void DoSave() {
    if (g_currentFilePath.empty()) {
        DoSaveAs(); // If no path, act like "Save As"
    } else {
        SCircuit ckt_data = createSaveDataFromCanvas(false);
        saveProject(g_currentFilePath, ckt_data);
    }
}

void DoOpenFile(const std::string& path) {
    SCircuit ckt_data;
    if (loadProject(path, ckt_data)) {
        loadCanvasFromSaveData(ckt_data);
        g_currentFilePath = path;
        recent_add(g_currentFilePath); // Add/move to top of recent list
    } else {
        std::cerr << "Failed to load project file: " << path << std::endl;
        // Optional: Add an error pop-up for the user here
    }
}

void DoOpen() {
    auto selection = pfd::open_file("Open Circuit", "", { "JSON Files", "*.json" }).result();
    if (!selection.empty()) {
        DoOpenFile(selection[0]);
    }
}

void RenderMenuBar() {
    if (ImGui::BeginMainMenuBar()) {

        if (ImGui::BeginMenu("File")) {
            if (ImGui::MenuItem("New", "Ctrl+N")) { NewProject(); }
            if (ImGui::MenuItem("Open", "Ctrl+O")) { DoOpen(); }
            if (ImGui::MenuItem("Save", "Ctrl+S")) { DoSave(); }
            if (ImGui::MenuItem("Save As...")) { DoSaveAs(); }
            ImGui::Separator();

            // --- THIS IS THE CORRECTED BLOCK ---
            // Use BeginMenu to create a submenu, not MenuItem
            if (ImGui::BeginMenu("Recent Files")) {
                std::vector<std::string> recent_files = recent_list();
                if (recent_files.empty()) {
                    ImGui::MenuItem("(Empty)", NULL, false, false);
                } else {
                    for (const auto& file_path : recent_files) {
                        if (ImGui::MenuItem(file_path.c_str())) {
                            DoOpenFile(file_path);
                        }
                    }
                }
                ImGui::EndMenu(); // This is now correctly paired with BeginMenu
            }

            ImGui::Separator();
            if (ImGui::MenuItem("Exit")) { exit(0); }
            ImGui::EndMenu();
        }

        // --- RUN & DEBUG MENUS (Unchanged) ---
        if (ImGui::BeginMenu("Run")) {
            if (ImGui::MenuItem("Run Simulation...")) {
                g_showRunPopup = true;
            }
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Debug")) {
            ImGui::MenuItem("Show Node Info", NULL, &g_showDebugInfo);
            ImGui::EndMenu();
        }

        RenderSubcircuitMenu();
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
    if (ImGui::Button("VPULSE")) { g_showVPulsePopup = true; } ImGui::SameLine();
    if (ImGui::Button("IPULSE")) { g_showIPulsePopup = true; } ImGui::SameLine();
    if (ImGui::Button("VDelta")) { g_showVDeltaPopup = true; } ImGui::SameLine();
    if (ImGui::Button("IDelta")) { g_showIDeltaPopup = true; } ImGui::SameLine();
    if (ImGui::Button("VPhase")) { g_showVPhasePopup = true; } ImGui::SameLine();
    if (ImGui::Button("VCVS")) {
        g_currentTool = VCVS;
        g_showDepSourcePopup = true;
        g_depPopupError.clear();
    } ImGui::SameLine();
    if (ImGui::Button("VCCS")) {
        g_currentTool = VCCS;
        g_showDepSourcePopup = true;
        g_depPopupError.clear();
    } ImGui::SameLine();
    if (ImGui::Button("Sub")) {
        g_currentTool = SUBCIRCUIT_TOOL;
        g_subcircuitSelectionStep = 0;
    } ImGui::SameLine();
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
        ImGui::RadioButton("Phase Sweep", (int*)&g_runChoice, RUN_PHASE);
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
        } else if (g_runChoice == RUN_PHASE) {
            ImGui::InputText("Base Freq (Hz)", g_phaseBaseFreqText, 64);
            ImGui::InputText("Start Phase (deg)", g_phaseStartText, 64);
            ImGui::InputText("Stop Phase (deg)", g_phaseStopText, 64);
            ImGui::InputText("Num Steps", g_phaseNumStepsText, 64);
            ImGui::InputText("Probe Node", g_phaseProbeText, 64);
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
                } else if (g_runChoice == RUN_PHASE) {
                    double base_freq = parseNumber(g_phaseBaseFreqText);
                    double phase_start_deg = parseNumber(g_phaseStartText);
                    double phase_stop_deg = parseNumber(g_phaseStopText);
                    int num_steps = std::stoi(g_phaseNumStepsText);
                    g_phase_probe_node = g_phaseProbeText;

                    g_phase_angles.clear();
                    g_phase_mags.clear();

                    for (int i = 0; i <= num_steps; ++i) {
                        double phase_deg = phase_start_deg + (phase_stop_deg - phase_start_deg) * i / num_steps;
                        double phase_rad = phase_deg * M_PI / 180.0;

                        // Find the first VPHASE/VSIN source and set its phase
                        for (auto *e: g_circuit.elements) {
                            if (auto vs = dynamic_cast<VoltageSource *>(e)) {
                                vs->phase = phase_rad;
                                break; // Only modify the first source found
                            }
                        }

                        std::map<std::string, cd> V;
                        solveACAtFrequency(g_circuit, base_freq, V);
                        cd v_out = V.count(g_phase_probe_node) ? V.at(g_phase_probe_node) : 0;

                        g_phase_angles.push_back(phase_deg);
                        g_phase_mags.push_back(std::abs(v_out));
                    }
                    g_showPhasePlot = true;
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


// Interpolates to find the y-value on a signal for a given x-value (time)
double interpolate(const std::vector<double>& x_data, const std::vector<double>& y_data, double x) {
    if (x <= x_data.front()) return y_data.front();
    if (x >= x_data.back()) return y_data.back();

    auto it = std::lower_bound(x_data.begin(), x_data.end(), x);
    size_t i = std::distance(x_data.begin(), it);

    double x1 = x_data[i-1], y1 = y_data[i-1];
    double x2 = x_data[i], y2 = y_data[i];

    return y1 + (y2 - y1) * (x - x1) / (x2 - x1);
}


void RenderTransientPlotWindow() {
    if (!g_showTransientPlot) {
        g_showAddTrace = false;
        g_showSubtractTrace = false;
        return;
    }

    ImGui::SetNextWindowSize(ImVec2(700, 500), ImGuiCond_FirstUseEver);
    ImGui::Begin("Transient Analysis Results", &g_showTransientPlot);

    // --- 1. Prepare signals and math (unchanged) ---
    g_plot_signal_names.clear();
    g_plot_signal_map.clear();
    int current_idx = 0;
    for (const auto& pair : g_waves.V) {
        g_plot_signal_names.push_back(pair.first.c_str());
        g_plot_signal_map[current_idx++] = {pair.first, 'V'};
    }
    for (const auto& pair : g_waves.I) {
        g_plot_signal_names.push_back(pair.first.c_str());
        g_plot_signal_map[current_idx++] = {pair.first, 'I'};
    }
    if ((g_showAddTrace || g_showSubtractTrace) && g_mathSignal1_idx != -1 && g_mathSignal2_idx != -1) {
        auto signal1_info = g_plot_signal_map[g_mathSignal1_idx];
        auto signal2_info = g_plot_signal_map[g_mathSignal2_idx];
        const auto& data1 = (signal1_info.second == 'V') ? g_waves.V.at(signal1_info.first) : g_waves.I.at(signal1_info.first);
        const auto& data2 = (signal2_info.second == 'V') ? g_waves.V.at(signal2_info.first) : g_waves.I.at(signal2_info.first);

        if (data1.size() == data2.size()) {
            g_addTrace.resize(data1.size());
            g_subtractTrace.resize(data1.size());
            for (size_t i = 0; i < data1.size(); ++i) {
                g_addTrace[i] = data1[i] + data2[i];
                g_subtractTrace[i] = data1[i] - data2[i];
            }
        }
    }

    // --- 2. Math Channel Controls (unchanged) ---
    if (ImGui::CollapsingHeader("Math Channel")) {
        ImGui::Combo("Signal A", &g_mathSignal1_idx, g_plot_signal_names.data(), g_plot_signal_names.size());
        ImGui::Combo("Signal B", &g_mathSignal2_idx, g_plot_signal_names.data(), g_plot_signal_names.size());
        ImGui::Checkbox("Plot Sum (A+B)", &g_showAddTrace);
        ImGui::SameLine();
        ImGui::Checkbox("Plot Difference (A-B)", &g_showSubtractTrace);
    }

    if (ImPlot::BeginPlot("Waveforms")) {
        ImPlot::SetupAxes("Time (s)", "Value");

        // Plot original and calculated traces (unchanged)
        for (const auto& pair : g_waves.V) ImPlot::PlotLine(pair.first.c_str(), g_waves.t.data(), pair.second.data(), g_waves.t.size());
        for (const auto& pair : g_waves.I) ImPlot::PlotLine(pair.first.c_str(), g_waves.t.data(), pair.second.data(), g_waves.t.size());

        if (g_showAddTrace && !g_addTrace.empty()) {
            std::string label = std::string(g_plot_signal_names[g_mathSignal1_idx]) + "+" + std::string(g_plot_signal_names[g_mathSignal2_idx]);
            ImPlot::PlotLine(label.c_str(), g_waves.t.data(), g_addTrace.data(), g_waves.t.size());
        }
        if (g_showSubtractTrace && !g_subtractTrace.empty()) {
            std::string label = std::string(g_plot_signal_names[g_mathSignal1_idx]) + "-" + std::string(g_plot_signal_names[g_mathSignal2_idx]);
            ImPlot::PlotLine(label.c_str(), g_waves.t.data(), g_subtractTrace.data(), g_waves.t.size());
        }

        // --- B. Handle all hover and click interactions ---
        if (ImPlot::IsPlotHovered()) {
            // --- Accurate Hover Tooltip (now uses pixels) ---
            ImVec2 mouse_pixels = ImGui::GetMousePos();
            int closest_series = -1;
            int closest_point_idx = -1;
            float min_dist_sq = 100.0f; // 10 pixel radius

            int series_id = 0;
            auto find_closest_px = [&](const std::vector<double>& data) {
                for (int i = 0; i < data.size(); ++i) {
                    ImVec2 p_pixels = ImPlot::PlotToPixels({ g_waves.t[i], data[i] });
                    float dx = mouse_pixels.x - p_pixels.x;
                    float dy = mouse_pixels.y - p_pixels.y;
                    float dist_sq = dx*dx + dy*dy;
                    if (dist_sq < min_dist_sq) {
                        min_dist_sq = dist_sq;
                        closest_series = series_id;
                        closest_point_idx = i;
                    }
                }
                series_id++;
            };
            for(const auto& pair : g_waves.V) find_closest_px(pair.second);
            for(const auto& pair : g_waves.I) find_closest_px(pair.second);

            if (closest_point_idx != -1) {
                auto signal_info = g_plot_signal_map.at(closest_series);
                const auto &data = (signal_info.second == 'V') ? g_waves.V.at(signal_info.first) : g_waves.I.at(
                        signal_info.first);
                double found_x = g_waves.t.at(closest_point_idx);
                double found_y = data.at(closest_point_idx);

                ImPlot::GetPlotDrawList()->AddCircleFilled(ImPlot::PlotToPixels({found_x, found_y}), 5,
                                                           IM_COL32(255, 0, 0, 255));
                ImGui::BeginTooltip();
                ImGui::Text("%s", signal_info.first.c_str());
                ImGui::Text("Time:  %.4f s", found_x);
                ImGui::Text("Value: %.4f", found_y);
                ImGui::EndTooltip();

                // Cursor Placement Logic
                if (ImGui::GetIO().KeyCtrl && ImGui::IsMouseClicked(ImGuiMouseButton_Middle)) {
                    g_cursor1_active = true;
                    g_cursor1_series_idx = closest_series;
                    g_cursor1_pos = { found_x, found_y };
                }
                if (ImGui::GetIO().KeyShift && ImGui::IsMouseClicked(ImGuiMouseButton_Middle)) {
                    g_cursor2_active = true;
                    g_cursor2_series_idx = closest_series;
                    g_cursor2_pos = { found_x, found_y };
                }
            }
        }

        // --- C. Draw and Constrain Draggable Cursors ---
        if (g_cursor1_active) {
            if (ImPlot::DragPoint(0, &g_cursor1_pos.x, &g_cursor1_pos.y, ImVec4(1,1,0,1), 4)) {
                // Point was dragged, now constrain it to the waveform
                auto info = g_plot_signal_map[g_cursor1_series_idx];
                const auto& data = (info.second == 'V') ? g_waves.V.at(info.first) : g_waves.I.at(info.first);
                g_cursor1_pos.y = interpolate(g_waves.t, data, g_cursor1_pos.x);
            }
        }
        if (g_cursor2_active) {
            if (ImPlot::DragPoint(1, &g_cursor2_pos.x, &g_cursor2_pos.y, ImVec4(0,1,1,1), 4)) {
                // Point was dragged, now constrain it to the waveform
                auto info = g_plot_signal_map[g_cursor2_series_idx];
                const auto& data = (info.second == 'V') ? g_waves.V.at(info.first) : g_waves.I.at(info.first);
                g_cursor2_pos.y = interpolate(g_waves.t, data, g_cursor2_pos.x);
            }
        }

        ImPlot::EndPlot();
    }
    ImGui::End();
}


void RenderACPlotWindow() {
    if (!g_showACPlot) return;

    ImGui::SetNextWindowSize(ImVec2(700, 500), ImGuiCond_FirstUseEver);
    ImGui::Begin("AC Sweep Results", &g_showACPlot);

    if (ImPlot::BeginPlot(("V(" + g_ac_probe_node + ")").c_str())) {
        // Setup axes with a log scale for Frequency
        ImPlot::SetupAxes("Frequency (Hz)", "Magnitude (dB)");
        ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);

        // Plot the main magnitude data
        if (!g_ac_freqs.empty()) {
            ImPlot::PlotLine("Magnitude", g_ac_freqs.data(), g_ac_mags_db.data(), g_ac_freqs.size());
        }

        // --- NEW: Cursor and Tooltip Logic ---
        if (ImPlot::IsPlotHovered()) {
            // Find the closest point using screen pixels for accuracy on the log scale
            ImVec2 mouse_pixels = ImGui::GetMousePos();
            int closest_point_idx = -1;
            float min_dist_sq = 100.0f; // 10 pixel radius

            if (!g_ac_freqs.empty()) {
                for (int i = 0; i < g_ac_freqs.size(); ++i) {
                    ImVec2 p_pixels = ImPlot::PlotToPixels({ g_ac_freqs[i], g_ac_mags_db[i] });
                    float dx = mouse_pixels.x - p_pixels.x;
                    float dy = mouse_pixels.y - p_pixels.y;
                    float dist_sq = dx*dx + dy*dy;
                    if (dist_sq < min_dist_sq) {
                        min_dist_sq = dist_sq;
                        closest_point_idx = i;
                    }
                }
            }

            // Show tooltip if a point is found
            if (closest_point_idx != -1) {
                double found_x = g_ac_freqs[closest_point_idx];
                double found_y = g_ac_mags_db[closest_point_idx];
                ImPlot::GetPlotDrawList()->AddCircleFilled(ImPlot::PlotToPixels({found_x, found_y}), 5, IM_COL32(255, 0, 0, 255));

                ImGui::BeginTooltip();
                ImGui::Text("Freq: %.2f Hz", found_x);
                ImGui::Text("Mag:  %.2f dB", found_y);
                ImGui::EndTooltip();

                // Place Cursors with Ctrl/Shift + Middle-Click
                if ((ImGui::GetIO().KeyCtrl || ImGui::GetIO().KeyShift) && ImGui::IsMouseClicked(ImGuiMouseButton_Middle)) {
                    ImPlotPoint found_point = { found_x, found_y };
                    if (ImGui::GetIO().KeyCtrl) { g_cursor1_active = true; g_cursor1_pos = found_point; }
                    if (ImGui::GetIO().KeyShift) { g_cursor2_active = true; g_cursor2_pos = found_point; }
                }
            }
        }

        // Draw and constrain draggable cursors
        if (g_cursor1_active) {
            if (ImPlot::DragPoint(0, &g_cursor1_pos.x, &g_cursor1_pos.y, ImVec4(1,1,0,1), 4)) {
                g_cursor1_pos.y = interpolate(g_ac_freqs, g_ac_mags_db, g_cursor1_pos.x);
            }
        }
        if (g_cursor2_active) {
            if (ImPlot::DragPoint(1, &g_cursor2_pos.x, &g_cursor2_pos.y, ImVec4(0,1,1,1), 4)) {
                g_cursor2_pos.y = interpolate(g_ac_freqs, g_ac_mags_db, g_cursor2_pos.x);
            }
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


void RenderVPulsePopup() {
    if (g_showVPulsePopup) {
        ImGui::OpenPopup("Pulse Voltage Source");
        g_showVPulsePopup = false;
    }

    if (ImGui::BeginPopupModal("Pulse Voltage Source", NULL, ImGuiWindowFlags_AlwaysAutoResize)) {
        ImGui::Text("Pulse Source Properties");
        ImGui::Separator();
        ImGui::InputText("Name", g_vpulseNameBuffer, 64);
        ImGui::InputText("Initial Value (V)", g_vpulseVInitial, 64);
        ImGui::InputText("Pulsed Value (V)", g_vpulseVOn, 64);
        ImGui::InputText("Delay Time (s)", g_vpulseTDelay, 64);
        ImGui::InputText("Rise Time (s)", g_vpulseTRise, 64);
        ImGui::InputText("Fall Time (s)", g_vpulseTFall, 64);
        ImGui::InputText("Pulse Width (s)", g_vpulseTOn, 64);
        ImGui::InputText("Period (s)", g_vpulseTPeriod, 64);

        if (ImGui::Button("OK", ImVec2(120, 0))) {
            g_currentTool = VPULSE;
            ImGui::CloseCurrentPopup();
        }
        ImGui::SameLine();
        if (ImGui::Button("Cancel", ImVec2(120, 0))) { ImGui::CloseCurrentPopup(); }
        ImGui::EndPopup();
    }
}


void RenderIPulsePopup() {
    if (g_showIPulsePopup) {
        ImGui::OpenPopup("Pulse Current Source");
        g_showIPulsePopup = false;
    }

    if (ImGui::BeginPopupModal("Pulse Current Source", NULL, ImGuiWindowFlags_AlwaysAutoResize)) {
        ImGui::Text("Pulse Source Properties");
        ImGui::Separator();
        ImGui::InputText("Name", g_ipulseNameBuffer, 64);
        ImGui::InputText("Initial Value (A)", g_ipulseIInitial, 64);
        ImGui::InputText("Pulsed Value (A)", g_ipulseIOn, 64);
        ImGui::InputText("Delay Time (s)", g_ipulseTDelay, 64);
        ImGui::InputText("Rise Time (s)", g_ipulseTRise, 64);
        ImGui::InputText("Fall Time (s)", g_ipulseTFall, 64);
        ImGui::InputText("Pulse Width (s)", g_ipulseTOn, 64);
        ImGui::InputText("Period (s)", g_ipulseTPeriod, 64);

        if (ImGui::Button("OK", ImVec2(120, 0))) {
            g_currentTool = IPULSE;
            ImGui::CloseCurrentPopup();
        }
        ImGui::SameLine();
        if (ImGui::Button("Cancel", ImVec2(120, 0))) { ImGui::CloseCurrentPopup(); }
        ImGui::EndPopup();
    }
}


void RenderVDeltaPopup() {
    if (g_showVDeltaPopup) {
        ImGui::OpenPopup("Delta Voltage Source");
        g_showVDeltaPopup = false;
    }
    if (ImGui::BeginPopupModal("Delta Voltage Source", NULL, ImGuiWindowFlags_AlwaysAutoResize)) {
        ImGui::Text("Delta Source Properties");
        ImGui::Separator();
        ImGui::InputText("Name", g_vdeltaNameBuffer, 64);
        ImGui::InputText("Pulse Time (s)", g_vdeltaTPulse, 64);
        ImGui::InputText("Area", g_vdeltaArea, 64);

        if (ImGui::Button("OK", ImVec2(120, 0))) {
            g_currentTool = VDELTA;
            ImGui::CloseCurrentPopup();
        }
        ImGui::SameLine();
        if (ImGui::Button("Cancel", ImVec2(120, 0))) { ImGui::CloseCurrentPopup(); }
        ImGui::EndPopup();
    }
}


void RenderIDeltaPopup() {
    if (g_showIDeltaPopup) {
        ImGui::OpenPopup("Delta Current Source");
        g_showIDeltaPopup = false;
    }
    if (ImGui::BeginPopupModal("Delta Current Source", NULL, ImGuiWindowFlags_AlwaysAutoResize)) {
        ImGui::Text("Delta Source Properties");
        ImGui::Separator();
        ImGui::InputText("Name", g_ideltaNameBuffer, 64);
        ImGui::InputText("Pulse Time (s)", g_ideltaTPulse, 64);
        ImGui::InputText("Area", g_ideltaArea, 64);

        if (ImGui::Button("OK", ImVec2(120, 0))) {
            g_currentTool = IDELTA;
            ImGui::CloseCurrentPopup();
        }
        ImGui::SameLine();
        if (ImGui::Button("Cancel", ImVec2(120, 0))) { ImGui::CloseCurrentPopup(); }
        ImGui::EndPopup();
    }
}


void RenderVPhasePopup() {
    if (g_showVPhasePopup) { ImGui::OpenPopup("Phase Voltage Source"); g_showVPhasePopup = false; }
    if (ImGui::BeginPopupModal("Phase Voltage Source", NULL, ImGuiWindowFlags_AlwaysAutoResize)) {
        ImGui::Text("Phase Source Properties");
        ImGui::Separator();
        ImGui::InputText("Name", g_vphaseNameBuffer, 64);
        ImGui::InputText("Amplitude (V)", g_vphaseAmpBuffer, 64);
        ImGui::InputText("Frequency (Hz)", g_vphaseFreqBuffer, 64);
        ImGui::InputText("Phase (deg)", g_vphasePhaseBuffer, 64);
        if (ImGui::Button("OK")) { g_currentTool = VPHASE; ImGui::CloseCurrentPopup(); }
        ImGui::SameLine();
        if (ImGui::Button("Cancel")) { ImGui::CloseCurrentPopup(); }
        ImGui::EndPopup();
    }
}


void RenderPhasePlotWindow() {
    if (!g_showPhasePlot) {
        g_showAddTrace = false;
        g_showSubtractTrace = false;
        return;
    }
    ImGui::SetNextWindowSize(ImVec2(700, 500), ImGuiCond_FirstUseEver);
    ImGui::Begin("Phase Sweep Results", &g_showPhasePlot);

    // --- 1. Prepare signals for Math Channel dropdowns ---
    // Note: This will only show one signal for now, but is ready for the future
    g_plot_signal_names.clear();
    g_plot_signal_map.clear();
    if (!g_phase_mags.empty()) {
        g_plot_signal_names.push_back(g_phase_probe_node.c_str());
        g_plot_signal_map[0] = {g_phase_probe_node, 'V'};
    }

    // --- 2. Math Channel Controls ---
    if (ImGui::CollapsingHeader("Math Channel")) {
        ImGui::Combo("Signal A", &g_mathSignal1_idx, g_plot_signal_names.data(), g_plot_signal_names.size());
        ImGui::Combo("Signal B", &g_mathSignal2_idx, g_plot_signal_names.data(), g_plot_signal_names.size());
        ImGui::Checkbox("Plot Sum (A+B)", &g_showAddTrace);
        ImGui::SameLine();
        ImGui::Checkbox("Plot Difference (A-B)", &g_showSubtractTrace);
    }

    // --- 3. Draw the Plot ---
    if (ImPlot::BeginPlot(("V(" + g_phase_probe_node + ")").c_str())) {
        ImPlot::SetupAxes("Phase (degrees)", "Magnitude (V)");

        // Plot the main signal
        if (!g_phase_angles.empty()) {
            ImPlot::PlotLine("Magnitude", g_phase_angles.data(), g_phase_mags.data(), g_phase_angles.size());
        }

        // --- Snap-to-Signal Cursor/Tooltip Logic ---
        if (ImPlot::IsPlotHovered()) {
            ImPlotPoint mouse = ImPlot::GetPlotMousePos();
            int closest_point_idx = -1;
            float min_dist_sq = 100.0f; // 10 pixel radius

            // Find the closest point on the magnitude curve
            if (!g_phase_mags.empty()) {
                ImVec2 mouse_pixels = ImGui::GetMousePos();
                for (int i = 0; i < g_phase_angles.size(); ++i) {
                    ImVec2 p_pixels = ImPlot::PlotToPixels({ g_phase_angles[i], g_phase_mags[i] });
                    float dx = mouse_pixels.x - p_pixels.x;
                    float dy = mouse_pixels.y - p_pixels.y;
                    float dist_sq = dx*dx + dy*dy;
                    if (dist_sq < min_dist_sq) {
                        min_dist_sq = dist_sq;
                        closest_point_idx = i;
                    }
                }
            }

            // Show tooltip
            if (closest_point_idx != -1) {
                double found_x = g_phase_angles[closest_point_idx];
                double found_y = g_phase_mags[closest_point_idx];
                ImPlot::GetPlotDrawList()->AddCircleFilled(ImPlot::PlotToPixels({found_x, found_y}), 5, IM_COL32(255, 0, 0, 255));
                ImGui::BeginTooltip();
                ImGui::Text("%s", g_phase_probe_node.c_str());
                ImGui::Text("Phase: %.2f deg", found_x);
                ImGui::Text("Magnitude: %.3f V", found_y);
                ImGui::EndTooltip();

                // Place Cursors
                if (ImGui::GetIO().KeyCtrl && ImGui::IsMouseClicked(ImGuiMouseButton_Middle)) {
                    g_cursor1_active = true;
                    g_cursor1_pos = { found_x, found_y };
                }
                if (ImGui::GetIO().KeyShift && ImGui::IsMouseClicked(ImGuiMouseButton_Middle)) {
                    g_cursor2_active = true;
                    g_cursor2_pos = { found_x, found_y };
                }
            }
        }

        // Draw and handle draggable cursors
        if (g_cursor1_active) {
            if (ImPlot::DragPoint(0, &g_cursor1_pos.x, &g_cursor1_pos.y, ImVec4(1,1,0,1), 4)) {
                g_cursor1_pos.y = interpolate(g_phase_angles, g_phase_mags, g_cursor1_pos.x);
            }
        }
        if (g_cursor2_active) {
            if (ImPlot::DragPoint(1, &g_cursor2_pos.x, &g_cursor2_pos.y, ImVec4(0,1,1,1), 4)) {
                g_cursor2_pos.y = interpolate(g_phase_angles, g_phase_mags, g_cursor2_pos.x);
            }
        }
        ImPlot::EndPlot();
    }
    ImGui::End();
}

void RenderDependentSourcePopup() {
    if (g_showDepSourcePopup) {
        ImGui::OpenPopup("Dependent Source Properties");
        g_showDepSourcePopup = false;
    }
    if (ImGui::BeginPopupModal("Dependent Source Properties", NULL, ImGuiWindowFlags_AlwaysAutoResize)) {
        ImGui::Text("Dependent Source Properties");
        ImGui::Separator();
        ImGui::InputText("Name", g_depNameBuffer, 64);
        ImGui::InputText("Gain", g_depGainBuffer, 64);
        ImGui::InputText("Controlling Node (+)", g_depCtrlPBuffer, 64);
        ImGui::InputText("Controlling Node (-)", g_depCtrlNBuffer, 64);

        // --- Display error message if it exists ---
        if (!g_depPopupError.empty()) {
            ImGui::TextColored(ImVec4(1.0f, 0.0f, 0.0f, 1.0f), "%s", g_depPopupError.c_str());
        }

        if (ImGui::Button("OK", ImVec2(120, 0))) {
            // --- NEW: Validation Logic ---
            std::string ctrlP = g_depCtrlPBuffer;
            std::string ctrlN = g_depCtrlNBuffer;

            if (ctrlP.empty() || ctrlN.empty()) {
                g_depPopupError = "Error: Control node names cannot be empty.";
            } else {
                g_depPopupError.clear();
                ImGui::CloseCurrentPopup(); // Success, close the popup
            }
        }
        ImGui::SameLine();
        if (ImGui::Button("Cancel", ImVec2(120, 0))) {
            g_currentTool = NONE;
            ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
    }
}


void RenderCursorWindow() {
    // Show this window only if at least one cursor is active
    if (!g_cursor1_active && !g_cursor2_active) {
        return;
    }

    // Begin a new window to display the measurements
    ImGui::Begin("Cursor Measurements", nullptr, ImGuiWindowFlags_AlwaysAutoResize);

    // --- Cursor 1 Data ---
    ImGui::Text("Cursor 1");
    ImGui::Separator();
    if (g_cursor1_active) {
        ImGui::Text("Time (Horz): %.4f s", g_cursor1_pos.x);
        ImGui::Text("Value (Vert): %.4f", g_cursor1_pos.y);
    } else {
        ImGui::Text("Time (Horz): --");
        ImGui::Text("Value (Vert): --");
    }

    ImGui::Spacing();

    // --- Cursor 2 Data ---
    ImGui::Text("Cursor 2");
    ImGui::Separator();
    if (g_cursor2_active) {
        ImGui::Text("Time (Horz): %.4f s", g_cursor2_pos.x);
        ImGui::Text("Value (Vert): %.4f", g_cursor2_pos.y);
    } else {
        ImGui::Text("Time (Horz): --");
        ImGui::Text("Value (Vert): --");
    }

    ImGui::Spacing();

    // --- Difference Data ---
    ImGui::Text("Difference");
    ImGui::Separator();
    if (g_cursor1_active && g_cursor2_active) {
        ImGui::Text("Horz (Δt): %.4f s", std::abs(g_cursor2_pos.x - g_cursor1_pos.x));
        ImGui::Text("Vert (ΔV): %.4f", std::abs(g_cursor2_pos.y - g_cursor1_pos.y));
    } else {
        ImGui::Text("Horz (Δt): --");
        ImGui::Text("Vert (ΔV): --");
    }

    ImGui::End();
}


void RenderSubcircuitCreationWindow() {
    if (g_subcircuitSelectionStep != 2) return;

    ImGui::OpenPopup("Create Subcircuit");
    if (ImGui::BeginPopupModal("Create Subcircuit", NULL, ImGuiWindowFlags_AlwaysAutoResize)) {
        ImGui::Text("External Node 1: %s", g_subcircuitNode1_name.c_str());
        ImGui::Text("External Node 2: %s", g_subcircuitNode2_name.c_str());
        ImGui::InputText("Subcircuit Name", g_subcircuitNameBuffer, 64);

        if (ImGui::Button("Create")) {
            SCircuit saved_data = createSaveDataFromCanvas(false); // Save without ground check
            g_subcircuitLibrary.push_back({
                                                  g_subcircuitNameBuffer,
                                                  saved_data.Snodes,
                                                  saved_data.Selements,
                                                  saved_data.Swires,
                                                  g_subcircuitNode1_name,
                                                  g_subcircuitNode2_name
                                          });

            SaveSubcircuitLibrary();

            std::cout << "Subcircuit '" << g_subcircuitNameBuffer << "' created." << std::endl;
            g_subcircuitSelectionStep = 0;
            g_currentTool = NONE;
            ImGui::CloseCurrentPopup();
        }
        ImGui::SameLine();
        if (ImGui::Button("Cancel")) {
            g_subcircuitSelectionStep = 0;
            g_currentTool = NONE;
            ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
    }
}

void RenderSubcircuitMenu() {
    // This adds a new "Sub" menu to the main menu bar
    if (ImGui::BeginMenu("Sub")) {
        if (g_subcircuitLibrary.empty()) {
            ImGui::MenuItem("(Empty)", NULL, false, false);
        } else {
            for (int i = 0; i < g_subcircuitLibrary.size(); ++i) {
                if (ImGui::MenuItem(g_subcircuitLibrary[i].name.c_str())) {
                    g_currentTool = SUBCIRCUIT_INSTANCE;
                    g_placingSubcircuit_idx = i;
                    g_previewIsVertical = false;
                }
            }
        }
        ImGui::EndMenu();
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
void DrawComponent(ImDrawList* drawList, const PlacedElement& element, int index, const ImVec2& canvas_p0) {
    // --- THIS IS THE MAIN CHANGE ---
    // If this element's index matches the selected index, draw it in blue
    const ImU32 color = (index == g_selectedIndex) ? IM_COL32(0, 100, 255, 255) : IM_COL32(0, 0, 0, 255);

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
            drawList->PathLineTo(ImVec2(pos.x - 7, pos.y - 7));
            drawList->PathBezierCubicCurveTo(ImVec2(pos.x - 3, pos.y - 11), ImVec2(pos.x + 3, pos.y - 3), ImVec2(pos.x + 7, pos.y - 7), 10);
            drawList->PathStroke(color, ImDrawFlags_None, thickness);
            break;
        }

        case VPULSE: {
            drawList->AddCircle(pos, 15.0f, color, 0, thickness);
            // Draw a small pulse/square wave symbol
            ImVec2 p = ImVec2(pos.x - 8, pos.y + 5);
            drawList->AddLine(p, ImVec2(p.x + 5, p.y), color, thickness);
            drawList->AddLine(ImVec2(p.x + 5, p.y), ImVec2(p.x + 5, p.y - 10), color, thickness);
            drawList->AddLine(ImVec2(p.x + 5, p.y - 10), ImVec2(p.x + 11, p.y - 10), color, thickness);
            drawList->AddLine(ImVec2(p.x + 11, p.y - 10), ImVec2(p.x + 11, p.y), color, thickness);
            drawList->AddLine(ImVec2(p.x + 11, p.y), ImVec2(p.x + 16, p.y), color, thickness);
            break;
        }

        case IPULSE: {
            drawList->AddCircle(pos, 15.0f, color, 0, thickness);
            if (element.isVertical) {
                drawList->AddLine(ImVec2(pos.x, pos.y - 10), ImVec2(pos.x, pos.y + 10), color, thickness);
                drawList->AddLine(ImVec2(pos.x, pos.y + 10), ImVec2(pos.x - 4, pos.y + 5), color, thickness);
                drawList->AddLine(ImVec2(pos.x, pos.y + 10), ImVec2(pos.x + 4, pos.y + 5), color, thickness);
            } else {
                drawList->AddLine(ImVec2(pos.x - 10, pos.y), ImVec2(pos.x + 10, pos.y), color, thickness);
                drawList->AddLine(ImVec2(pos.x + 10, pos.y), ImVec2(pos.x + 5, pos.y - 4), color, thickness);
                drawList->AddLine(ImVec2(pos.x + 10, pos.y), ImVec2(pos.x + 5, pos.y + 4), color, thickness);
            }
            ImVec2 p = ImVec2(pos.x - 7, pos.y - 7);
            drawList->AddLine(ImVec2(p.x, p.y+4), ImVec2(p.x, p.y), color, thickness);
            drawList->AddLine(ImVec2(p.x, p.y), ImVec2(p.x+7, p.y), color, thickness);
            drawList->AddLine(ImVec2(p.x+7, p.y), ImVec2(p.x+7, p.y+4), color, thickness);
            break;
        }

        case VDELTA: {
            drawList->AddCircle(pos, 15.0f, color, 0, thickness);
            // Draw a spike symbol for delta
            drawList->AddLine(ImVec2(pos.x, pos.y + 10), ImVec2(pos.x, pos.y - 2), color, thickness);
            drawList->AddLine(ImVec2(pos.x, pos.y - 2), ImVec2(pos.x - 4, pos.y + 3), color, thickness);
            drawList->AddLine(ImVec2(pos.x, pos.y - 2), ImVec2(pos.x + 4, pos.y + 3), color, thickness);
            break;
        }

        case IDELTA: {
            // Draw C-Source Body (Circle + Arrow)
            drawList->AddCircle(pos, 15.0f, color, 0, thickness);
            if (element.isVertical) {
                drawList->AddLine(ImVec2(pos.x, pos.y - 10), ImVec2(pos.x, pos.y + 10), color, thickness);
                drawList->AddLine(ImVec2(pos.x, pos.y + 10), ImVec2(pos.x - 4, pos.y + 5), color, thickness);
                drawList->AddLine(ImVec2(pos.x, pos.y + 10), ImVec2(pos.x + 4, pos.y + 5), color, thickness);
            } else {
                drawList->AddLine(ImVec2(pos.x - 10, pos.y), ImVec2(pos.x + 10, pos.y), color, thickness);
                drawList->AddLine(ImVec2(pos.x + 10, pos.y), ImVec2(pos.x + 5, pos.y - 4), color, thickness);
                drawList->AddLine(ImVec2(pos.x + 10, pos.y), ImVec2(pos.x + 5, pos.y + 4), color, thickness);
            }
            // Draw a spike symbol for delta
            drawList->AddLine(ImVec2(pos.x, pos.y - 7), ImVec2(pos.x, pos.y + 3), color, thickness);
            drawList->AddLine(ImVec2(pos.x, pos.y - 7), ImVec2(pos.x - 3, pos.y - 2), color, thickness);
            drawList->AddLine(ImVec2(pos.x, pos.y - 7), ImVec2(pos.x + 3, pos.y - 2), color, thickness);
            break;
        }

        case VPHASE: {
            drawList->AddCircle(pos, 15.0f, color, 0, thickness);
            drawList->PathLineTo(ImVec2(pos.x - 10, pos.y));
            drawList->PathBezierCubicCurveTo(ImVec2(pos.x - 5, pos.y - 10), ImVec2(pos.x + 5, pos.y + 10), ImVec2(pos.x + 10, pos.y), 10);
            drawList->PathStroke(color, ImDrawFlags_None, thickness);
            drawList->AddText(ImVec2(pos.x - 4, pos.y - 15), color, "f");
            break;
        }

        case VCVS: {
            ImVec2 p1_d = {pos.x - 15, pos.y}, p2_d = {pos.x, pos.y - 15};
            ImVec2 p3_d = {pos.x + 15, pos.y}, p4_d = {pos.x, pos.y + 15};
            if (element.isVertical) {
                p1_d = {pos.x, pos.y - 15}; p2_d = {pos.x + 15, pos.y};
                p3_d = {pos.x, pos.y + 15}; p4_d = {pos.x - 15, pos.y};
            }
            drawList->AddQuad(p1_d, p2_d, p3_d, p4_d, color, thickness);

            // Add '+' and '-' signs to show polarity
            if (element.isVertical) {
                drawList->AddLine({pos.x - 3, pos.y - 7}, {pos.x + 3, pos.y - 7}, color, thickness);
                drawList->AddLine({pos.x, pos.y - 9}, {pos.x, pos.y - 5}, color, thickness);
                drawList->AddLine({pos.x - 3, pos.y + 7}, {pos.x + 3, pos.y + 7}, color, thickness);
            } else {
                drawList->AddLine({pos.x - 9, pos.y}, {pos.x - 5, pos.y}, color, thickness);
                drawList->AddLine({pos.x - 7, pos.y - 2}, {pos.x - 7, pos.y + 2}, color, thickness);
                drawList->AddLine({pos.x + 5, pos.y}, {pos.x + 9, pos.y}, color, thickness);
            }
            break;
        }

        case VCCS: {
            // Diamond shape
            ImVec2 p1_d = {pos.x - 15, pos.y}, p2_d = {pos.x, pos.y - 15};
            ImVec2 p3_d = {pos.x + 15, pos.y}, p4_d = {pos.x, pos.y + 15};
            if (element.isVertical) {
                p1_d = {pos.x, pos.y - 15}; p2_d = {pos.x + 15, pos.y};
                p3_d = {pos.x, pos.y + 15}; p4_d = {pos.x - 15, pos.y};
            }
            drawList->AddQuad(p1_d, p2_d, p3_d, p4_d, color, thickness);

            // Arrow for Current source
            drawList->AddLine({pos.x - 5, pos.y}, {pos.x + 5, pos.y}, color, thickness);
            drawList->AddLine({pos.x + 5, pos.y}, {pos.x, pos.y - 4}, color, thickness);
            drawList->AddLine({pos.x + 5, pos.y}, {pos.x, pos.y + 4}, color, thickness);
            break;
        }

        case SUBCIRCUIT_INSTANCE: {
            float half_w = element.size.x / 2.0f;
            float half_h = element.size.y / 2.0f;
            if(element.isVertical) std::swap(half_w, half_h);
            drawList->AddRect({pos.x - half_w, pos.y - half_h}, {pos.x + half_w, pos.y + half_h}, color, 4.0f, 0, thickness);
            if (element.subcircuit_def_idx != -1) {
                const char* label = g_subcircuitLibrary[element.subcircuit_def_idx].name.c_str();
                ImVec2 text_size = ImGui::CalcTextSize(label);
                drawList->AddText({pos.x - text_size.x / 2, pos.y - text_size.y / 2}, color, label);
            }
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
        } else if (element.type == VPULSE) {
            std::string offset_str = formatValueWithSI(element.v_initial);
            std::string high_str = formatValueWithSI(element.v_on);
            std::string period_str = formatValueWithSI(element.t_period);
            std::string on_str = formatValueWithSI(element.t_on); // Using Ton as "dutyperiod"
            snprintf(display_str, 128, "VPULSE(%s %s %s %s)", offset_str.c_str(), high_str.c_str(), period_str.c_str(), on_str.c_str());
        } else if (element.type == IPULSE) {
            std::string offset_str = formatValueWithSI(element.i_initial);
            std::string high_str = formatValueWithSI(element.i_on);
            std::string period_str = formatValueWithSI(element.t_period);
            std::string on_str = formatValueWithSI(element.t_on);
            snprintf(display_str, 128, "IPULSE(%s %s %s %s)", offset_str.c_str(), high_str.c_str(), period_str.c_str(), on_str.c_str());
        } else if (element.type == VCVS || element.type == VCCS) {
            snprintf(display_str, 128, "(%s %s %.2f)", element.ctrlNodeP_name, element.ctrlNodeN_name, element.gain);
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


void HandleCanvasShortcuts() {
    if (!ImGui::IsWindowHovered(ImGuiHoveredFlags_RootAndChildWindows) || ImGui::IsAnyItemActive()) {
        return;
    }

    // --- Action Keys (Highest Priority) ---
    if (ImGui::IsKeyPressed(ImGuiKey_Escape)) {
        if (g_isMovingElement) {
            g_placedElements[g_selectedIndex].pos = g_originalPos;
            g_isMovingElement = false;
            g_selectedIndex = -1;
        } else {
            g_currentTool = NONE;
            g_isPlacingWire = false;
            g_previewIsVertical = false;
        }
        return;
    }

    if (ImGui::IsKeyPressed(ImGuiKey_Delete) && g_selectedIndex != -1) {
        g_placedElements.erase(g_placedElements.begin() + g_selectedIndex);
        g_selectedIndex = -1;
        g_isMovingElement = false;
        return;
    }

    // --- Rotation (Second Priority) ---
    if (ImGui::GetIO().KeyCtrl && ImGui::IsKeyPressed(ImGuiKey_R)) {
        // Priority 1: Rotate a selected element
        if (g_selectedIndex != -1) {
            g_placedElements[g_selectedIndex].isVertical = !g_placedElements[g_selectedIndex].isVertical;
        }
            // Priority 2: Rotate the current tool's preview
        else if (g_currentTool != NONE && g_currentTool != WIRE) {
            g_previewIsVertical = !g_previewIsVertical;
        }
        return; // Action handled
    }

    // --- Tool Selection (Lowest Priority) ---
    // This part only runs if the keys above weren't pressed.
    auto request_component_popup = [](ToolType type) {
        g_pendingTool = type; g_showComponentPopup = true;
        g_componentNameBuffer[0] = '\0'; g_componentValueBuffer[0] = '\0';
        g_previewIsVertical = false; g_componentPopupError.clear();
    };

    if (ImGui::IsKeyPressed(ImGuiKey_R)) { request_component_popup(RESISTOR); }
    else if (ImGui::IsKeyPressed(ImGuiKey_C)) { request_component_popup(CAPACITOR); }
    else if (ImGui::IsKeyPressed(ImGuiKey_L)) { request_component_popup(INDUCTOR); }
    else if (ImGui::IsKeyPressed(ImGuiKey_V)) { request_component_popup(VSOURCE); }
    else if (ImGui::IsKeyPressed(ImGuiKey_I)) { request_component_popup(CSOURCE); }
    else if (ImGui::IsKeyPressed(ImGuiKey_D)) { request_component_popup(DIODE); }
    else if (ImGui::IsKeyPressed(ImGuiKey_G)) { g_currentTool = GROUND; g_previewIsVertical = false; }
    else if (ImGui::IsKeyPressed(ImGuiKey_W)) { g_currentTool = WIRE; g_isPlacingWire = false; }
}


bool IsPointOnWire(const ImVec2& p, const Wire& wire) {
    const float threshold = 3.0f; // How close the click must be
    float dx = wire.p2.x - wire.p1.x;
    float dy = wire.p2.y - wire.p1.y;

    // If the wire is just a point, we can't connect to its middle
    if (dx == 0 && dy == 0) return false;

    // Check if the point is within the bounding box of the wire segment
    if (p.x < std::min(wire.p1.x, wire.p2.x) - threshold || p.x > std::max(wire.p1.x, wire.p2.x) + threshold ||
        p.y < std::min(wire.p1.y, wire.p2.y) - threshold || p.y > std::max(wire.p1.y, wire.p2.y) + threshold) {
        return false;
    }

    // Calculate the distance from the point to the line
    float dist = std::abs(dx * (wire.p1.y - p.y) - (wire.p1.x - p.x) * dy) / std::sqrt(dx * dx + dy * dy);

    return dist < threshold;
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

    HandleCanvasShortcuts();

    if (ImGui::IsWindowHovered(ImGuiHoveredFlags_RootAndChildWindows)) {
        const ImVec2 mouse_pos_in_canvas = ImVec2(ImGui::GetMousePos().x - canvas_p0.x, ImGui::GetMousePos().y - canvas_p0.y);
        const ImVec2 snapped_pos = SnapToGrid(mouse_pos_in_canvas, GRID_STEP);

        // --- NEW: Logic to move the selected element ---
        if (g_isMovingElement && g_selectedIndex != -1) {
            g_placedElements[g_selectedIndex].pos = snapped_pos;
        }

        // --- Placement Logic ---
        if (ImGui::IsMouseClicked(ImGuiMouseButton_Left) || ImGui::IsKeyPressed(ImGuiKey_Enter)) {

            if (g_currentTool == SUBCIRCUIT_TOOL) {
                buildCircuit(false); // Analyze connections without requiring a ground
                int closest_connector_id = -1;
                float min_dist_sq = 100.0f;
                for (const auto& conn : g_allConnectors) {
                    float dx = conn.pos.x - snapped_pos.x; float dy = conn.pos.y - snapped_pos.y;
                    if ((dx * dx + dy * dy) < min_dist_sq) {
                        min_dist_sq = (dx * dx + dy * dy);
                        closest_connector_id = conn.id;
                    }
                }

                if (closest_connector_id != -1) {
                    std::string clicked_node_name = g_connectorToNodeName.at(closest_connector_id);
                    if (g_subcircuitSelectionStep == 0) {
                        g_subcircuitNode1_name = clicked_node_name;
                        g_subcircuitNode1_pos = snapped_pos;
                        g_subcircuitSelectionStep = 1;
                    } else if (g_subcircuitSelectionStep == 1 && clicked_node_name != g_subcircuitNode1_name) {
                        g_subcircuitNode2_name = clicked_node_name;
                        g_subcircuitNode2_pos = snapped_pos;
                        g_subcircuitSelectionStep = 2; // Both nodes selected, popup will now appear
                    }
                }
            }
            else if (g_currentTool == NONE) {
                if (g_isMovingElement) {
                    // If we were moving an element, this click places it.
                    g_isMovingElement = false;
                    g_selectedIndex = -1;
                } else {
                    // If no tool is active, try to select an element.
                    g_selectedIndex = -1; // Deselect first
                    for (int i = g_placedElements.size() - 1; i >= 0; --i) {
                        const auto& elem = g_placedElements[i];
                        ImVec2 min = ImVec2(elem.pos.x - PIN_OFFSET, elem.pos.y - 15);
                        ImVec2 max = ImVec2(elem.pos.x + PIN_OFFSET, elem.pos.y + 15);
                        if (elem.isVertical) {
                            min = ImVec2(elem.pos.x - 15, elem.pos.y - PIN_OFFSET);
                            max = ImVec2(elem.pos.x + 15, elem.pos.y + PIN_OFFSET);
                        }

                        ImVec2 rect_min_abs = ImVec2(canvas_p0.x + min.x, canvas_p0.y + min.y);
                        ImVec2 rect_max_abs = ImVec2(canvas_p0.x + max.x, canvas_p0.y + max.y);
                        if (ImGui::IsMouseHoveringRect(rect_min_abs, rect_max_abs)) {
                            g_selectedIndex = i;
                            g_isMovingElement = true;
                            g_originalPos = elem.pos;
                            break;
                        }
                    }
                }
            } else {
                if (g_currentTool == VOLTMETER) {
                    buildCircuit();
                    int closest_connector_id = -1;
                    float min_dist_sq = 100.0f;
                    for (const auto &conn: g_allConnectors) {
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
                        } catch (const std::exception &e) {
                            g_simulationError = e.what();
                            g_showDCResults = true;
                        }
                    }
                    g_currentTool = NONE;
                } else if (g_currentTool == WIRE) {
                    if (!g_isPlacingWire) {
                        // First click: Start the wire
                        g_isPlacingWire = true;
                        g_wireStartPos = snapped_pos;
                    } else {
                        // Second click: Finish the wire and check for junctions
                        ImVec2 end_pos = snapped_pos;

                        // Check if the end point lands on an existing wire
                        int split_wire_index = -1;
                        for (int i = 0; i < g_wires.size(); ++i) {
                            if (IsPointOnWire(end_pos, g_wires[i])) {
                                split_wire_index = i;
                                break;
                            }
                        }

                        if (split_wire_index != -1) {
                            // We found a T-junction. Split the existing wire.
                            Wire old_wire = g_wires[split_wire_index];

                            // Remove the old, long wire
                            g_wires.erase(g_wires.begin() + split_wire_index);

                            // Add two new, shorter wires in its place
                            g_wires.push_back({old_wire.p1, end_pos});
                            g_wires.push_back({end_pos, old_wire.p2});
                        }

                        // Add the new wire the user just drew
                        g_wires.push_back({g_wireStartPos, end_pos});

                        g_isPlacingWire = false;
                        g_currentTool = NONE; // Deselect the tool
                    }
                } else if (g_currentTool != NONE) {
                    PlacedElement new_elem{};
                    new_elem.type = g_currentTool;
                    new_elem.pos = snapped_pos;
                    new_elem.isVertical = g_previewIsVertical;

                    if (g_currentTool == SUBCIRCUIT_INSTANCE) {
                        new_elem.name = g_subcircuitLibrary.at(g_placingSubcircuit_idx).name;
                        new_elem.subcircuit_def_idx = g_placingSubcircuit_idx;
                        new_elem.size = ImVec2(80, 50); // Set the default size
                    } else if (g_currentTool == VSIN) {
                        new_elem.name = g_vacNameBuffer;
                        new_elem.value = parseNumber(g_vacDcOffsetText);
                        new_elem.amp = parseNumber(g_vacAmpText);
                        new_elem.freq = parseNumber(g_vacFreqText);
                    } else if (g_currentTool == CSIN) {
                        new_elem.name = g_iacNameBuffer;
                        new_elem.value = parseNumber(g_iacDcOffsetText);
                        new_elem.amp = parseNumber(g_iacAmpText);
                        new_elem.freq = parseNumber(g_iacFreqText);
                    } else if (g_currentTool == VPULSE) {
                        new_elem.name = g_vpulseNameBuffer;
                        new_elem.v_initial = parseNumber(g_vpulseVInitial);
                        new_elem.v_on = parseNumber(g_vpulseVOn);
                        new_elem.t_delay = parseNumber(g_vpulseTDelay);
                        new_elem.t_rise = parseNumber(g_vpulseTRise);
                        new_elem.t_fall = parseNumber(g_vpulseTFall);
                        new_elem.t_on = parseNumber(g_vpulseTOn);
                        new_elem.t_period = parseNumber(g_vpulseTPeriod);
                    } else if (g_currentTool == IPULSE) {
                        new_elem.name = g_ipulseNameBuffer;
                        new_elem.i_initial = parseNumber(g_ipulseIInitial);
                        new_elem.i_on = parseNumber(g_ipulseIOn);
                        new_elem.t_delay = parseNumber(g_ipulseTDelay);
                        new_elem.t_rise = parseNumber(g_ipulseTRise);
                        new_elem.t_fall = parseNumber(g_ipulseTFall);
                        new_elem.t_on = parseNumber(g_ipulseTOn);
                        new_elem.t_period = parseNumber(g_ipulseTPeriod);
                    } else if (g_currentTool == VDELTA) {
                        new_elem.name = g_vdeltaNameBuffer;
                        new_elem.t_pulse = parseNumber(g_vdeltaTPulse);
                        new_elem.area = parseNumber(g_vdeltaArea);
                    } else if (g_currentTool == IDELTA) {
                        new_elem.name = g_ideltaNameBuffer;
                        new_elem.t_pulse = parseNumber(g_ideltaTPulse);
                        new_elem.area = parseNumber(g_ideltaArea);
                    } else if (g_currentTool == VPHASE) {
                        new_elem.name = g_vphaseNameBuffer;
                        new_elem.amp = parseNumber(g_vphaseAmpBuffer);
                        new_elem.freq = parseNumber(g_vphaseFreqBuffer);
                        new_elem.phase = parseNumber(g_vphasePhaseBuffer) * M_PI / 180.0; // Convert to radians
                    } else if (g_currentTool == VCVS || g_currentTool == VCCS) {
                        new_elem.name = g_depNameBuffer;
                        new_elem.gain = parseNumber(g_depGainBuffer);
                        strcpy_s(new_elem.ctrlNodeP_name, g_depCtrlPBuffer);
                        strcpy_s(new_elem.ctrlNodeN_name, g_depCtrlNBuffer);
                    } else { // Regular components (R, C, L, D, I)
                        new_elem.name = g_componentNameBuffer;
                        new_elem.value = parseNumber(g_componentValueBuffer);
                    }

                    if (g_currentTool == GROUND) {
                        new_elem.name = "GND";
                        new_elem.value = 0.0;
                    }

                    g_placedElements.push_back(new_elem);

                    g_currentTool = NONE;
                    g_previewIsVertical = false;
                }
            }
        }


        // --- Preview Logic ---
        // --- Preview Logic ---
        if (g_currentTool != NONE) {
            if (g_currentTool == WIRE && g_isPlacingWire) {
                ImVec2 start_abs = ImVec2(canvas_p0.x + g_wireStartPos.x, canvas_p0.y + g_wireStartPos.y);
                ImVec2 end_abs = ImVec2(canvas_p0.x + snapped_pos.x, canvas_p0.y + snapped_pos.y);
                drawList->AddLine(start_abs, end_abs, IM_COL32(0, 0, 255, 255), 1.5f);
            }
            else if (g_currentTool == VOLTMETER) {
                drawList->AddCircle(ImVec2(canvas_p0.x + snapped_pos.x, canvas_p0.y + snapped_pos.y), 15.0f, IM_COL32(23, 107, 135, 255), 0, 1.5f);
                drawList->AddText(ImVec2(canvas_p0.x + snapped_pos.x - 4, canvas_p0.y + snapped_pos.y - 8), IM_COL32(23, 107, 135, 255), "V");
            }
                // --- NEW: Special preview for Subcircuit Instances ---
            else if (g_currentTool == SUBCIRCUIT_INSTANCE) {
                PlacedElement preview_element{};
                preview_element.type = SUBCIRCUIT_INSTANCE;
                preview_element.pos = snapped_pos;
                preview_element.isVertical = g_previewIsVertical;
                preview_element.subcircuit_def_idx = g_placingSubcircuit_idx; // Pass the index
                preview_element.size = ImVec2(80, 50); // Use the default size
                DrawComponent(drawList, preview_element, -1, canvas_p0);
            }
                // This condition now correctly excludes all special tools
            else if (g_currentTool != SUBCIRCUIT_TOOL) {
                PlacedElement preview_element = { g_currentTool, snapped_pos, "", 0.0, g_previewIsVertical };
                DrawComponent(drawList, preview_element, -1, canvas_p0);
            }
        }
    }

    // --- Drawing Placed Items ---
    // FIX: This is the single, correct loop. The second one has been removed.
    for (int i = 0; i < g_placedElements.size(); ++i) {
        DrawComponent(drawList, g_placedElements[i], i, canvas_p0);
    }
    for (const auto& wire : g_wires) {
        ImVec2 p1 = ImVec2(canvas_p0.x + wire.p1.x, canvas_p0.y + wire.p1.y);
        ImVec2 p2 = ImVec2(canvas_p0.x + wire.p2.x, canvas_p0.y + wire.p2.y);
        drawList->AddLine(p1, p2, IM_COL32(0, 0, 0, 255), 1.5f);
        drawList->AddRectFilled(ImVec2(p1.x - 2, p1.y - 2), ImVec2(p1.x + 2, p1.y + 2), IM_COL32(0, 0, 0, 255));
        drawList->AddRectFilled(ImVec2(p2.x - 2, p2.y - 2), ImVec2(p2.x + 2, p2.y + 2), IM_COL32(0, 0, 0, 255));
    }
    for (const auto& probe : g_probes) DrawProbe(drawList, probe, canvas_p0);


    if (g_currentTool == SUBCIRCUIT_TOOL) {
        const char* instruction = "Select the FIRST external node (+ Port).";
        if (g_subcircuitSelectionStep == 1) {
            instruction = "Select the SECOND external node (- Port).";
            // Draw a marker on the first selected port
            drawList->AddText(ImVec2(canvas_p0.x + g_subcircuitNode1_pos.x, canvas_p0.y + g_subcircuitNode1_pos.y - 20), IM_COL32(0, 150, 0, 255), "+ Port");
        }
        // Draw the main instruction text
        drawList->AddText(ImVec2(canvas_p0.x + 10, canvas_p0.y + 10), IM_COL32(255, 0, 0, 255), instruction);
    }


    UpdateAndDrawDebugInfo(drawList, canvas_p0);

    ImGui::EndChild();
}


void HandleGlobalShortcuts() {
    // Don't process shortcuts if a text input is active
    if (ImGui::IsAnyItemActive()) {
        return;
    }

    const ImGuiIO& io = ImGui::GetIO();
    if (io.KeyCtrl) {
        if (ImGui::IsKeyPressed(ImGuiKey_N)) { NewProject(); }
        if (ImGui::IsKeyPressed(ImGuiKey_O)) { DoOpen(); }
        if (ImGui::IsKeyPressed(ImGuiKey_S)) { DoSave(); }
    }
}

// *********************************************************************************************************************************
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

    GLFWwindow* window = glfwCreateWindow(1280, 720, "New Gen z Circuit Simulator(V2.2.1)", NULL, NULL);
    if (!window) {
        std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);

    {
        int width, height, channels;
        // Make sure you have an "icon.png" file in your project folder
        unsigned char* pixels = stbi_load("icon.png", &width, &height, &channels, 4);
        if (pixels) {
            GLFWimage images[1];
            images[0].width = width;
            images[0].height = height;
            images[0].pixels = pixels;
            glfwSetWindowIcon(window, 1, images);
            stbi_image_free(pixels);
        } else {
            std::cerr << "Error: Could not load icon.png" << std::endl;
        }
    }

    glfwSwapInterval(1); // Enable vsync

    // 2. Initialize GLEW
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW" << std::endl;
        return -1;
    }

    // 3. Setup Dear ImGui and ImPlot contexts
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    LoadSubcircuitLibrary();
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
        HandleGlobalShortcuts();

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
        RenderVPulsePopup();
        RenderIPulsePopup();
        RenderVDeltaPopup();
        RenderIDeltaPopup();
        RenderVPhasePopup();
        RenderDependentSourcePopup();
        RenderCursorWindow();
        RenderSubcircuitCreationWindow();
        RenderRunPopup();
        RenderDCResultsWindow();
        RenderTransientPlotWindow();
        RenderACPlotWindow();
        RenderPhasePlotWindow();
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