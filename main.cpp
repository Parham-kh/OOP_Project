#include <bits/stdc++.h>
#include <unordered_map>
#include <filesystem>
#include <regex>


using namespace std;
//---------------------------------------------------------------------------------------------------------------------------------------
class Error : public exception {
    string message;

public:
    Error(const string& msg) : message(msg) {}
    const char* what() const noexcept override {
        return message.c_str();
    }
};


//---------------------------------------------------------------------------------------------------------------------------------------
bool success = false;

namespace Model {
    class Node {
    private:
        string name="";
        double voltage=0;

    public:
        Node(string name, double voltage) :name(name),voltage(voltage){}

        string getName() { return name; }

        double getVoltage() {
            return voltage;
        }
    };


    //---------------------------------------------------------------------------------------------------------------------------------------
    class Element {
    private:
        Node* node1;
        Node* node2;

    protected:
        string name="";
        string type="";
        double voltage=0;
        double current=0;

    public:
        Element(string name, string type, Node* node1, Node* node2, double voltage,double current)
                : name(name), type(type), node1(node1), node2(node2), voltage(voltage), current(current) {}

        string getName() const {
            return name;
        }

        string getType() const {
            return type;
        }

        Node* getNode1() {
            return node1;
        }

        Node* getNode2() {
            return node2;
        }

        virtual double getValue() {}

        virtual ~Element() = default;
    };


    class Resistor:public Element {
    private:
        double resistance = 0;

    public:
        Resistor(string name, Node* node1, Node* node2, double resistance)
                : Element(name, "Resistor", node1, node2, 0, 0), resistance(resistance) {}

        double getValue() override {
            return resistance;
        }
    };


    class Inductor:public Element{
    private:
        double inductance = 0;

    public:
        Inductor(string name, Node* node1, Node* node2, double inductance)
                : Element(name, "Inductor", node1, node2, 0, 0), inductance(inductance) {}

        double getValue() override {
            return inductance;
        }
    };


    class Capacitor : public Element {
    private:
        double capacitance = 0;

    public:
        Capacitor(string name, Node* node1, Node* node2, double capacitance)
                : Element(name, "Capacitor", node1, node2, 0, 0), capacitance(capacitance) {}

        double getValue() override {
            return capacitance;
        }
    };


    class CurrentSource:public Element{
    private:
        double current = 1;

    public:
        CurrentSource(string name, Node* node1, Node* node2, double current)
                : Element(name, "CurrentSource", node1, node2, 0, 0), current(current) {}

        double getValue() override {
            return current;
        }
    };


    class VoltageSource: public Element{
    private:
        double voltage = 1;

    public:
        VoltageSource(string name, Node* node1, Node* node2, double voltage)
                : Element(name, "VoltageSource", node1, node2, 0, 0), voltage(voltage) {}

        double getValue() override {
            return voltage;
        }
    };
};


//---------------------------------------------------------------------------------------------------------------------------------------
using namespace Model;
namespace Controller {
    struct TimePointResult {
        double time;
        unordered_map<string, double> nodeVoltages;
        unordered_map<string, double> elementCurrents;
    };

    class Circuit {
    private:
        vector<Element*> elements;
        vector<Node*> nodes;
        vector<TimePointResult> simulationResults;

    public:
        vector<Node*> getNodes() {
            return nodes;
        }


        double Coefficient_management(string a) {
            regex r1("(\\S*)([k|m|u|n|Meg])");
            smatch match;
            if(regex_search(a, match, r1)) {
                double output = stod(match[1].str());
                if(match[2]=="k") {
                    output *=1000;
                }
                else if(match[2]=="m") {
                    output *=0.001;
                }
                else if(match[2]=="u") {
                    output *=0.000001;
                }
                else if(match[2]=="n") {
                    output *=0.000000001;
                }
                else if(match[2]=="Meg") {
                    output *=1000000;
                }
                return output;
            }
            else {
                return stod(a);
            }
        }

        void addResistor(string name,string n1,string n2,string resistance) {
            try {
                if(name[0]!='R') {
                    throw Error("Error: Element "+name+" not found in library");
                }
                for(auto e:elements) {
                    if(e->getName()==name) {
                        throw Error("Error: Resistor "+ name +" already exists in the circuit");
                    }
                }
                Node* node1=new Node(n1,0);
                nodes.push_back(node1);

                Node* node2=new Node(n2,0);
                nodes.push_back(node2);
                Resistor* resistor=new Resistor(name,node1,node2,Coefficient_management(resistance));
                elements.push_back(resistor);
                success = true;
            } catch(const exception& e) {
                cout<<e.what()<<endl;
            }
        }

        void deleteResistor(string name) {
            bool found=false;
            for(int i=0;i<elements.size();i++) {
                if(elements[i]->getName() == name) {
                    elements.erase(elements.begin()+i);
                    found=true;
                }
            }
            if(!found) {
                throw Error("Error: Cannot delete resistor; component not found");
            }
        }

        void addCapacitor(string name,string n1,string n2,string capacitance) {
            try {
                for(auto e : elements) {
                    if(e->getName() == name) {
                        throw Error("Error: Capacitor "+ name +" already exists in the circuit");
                    }
                }
                Node* node1 = new Node(n1,0);
                Node* node2 = new Node(n2,0);
                nodes.push_back(node1);
                nodes.push_back(node2);
                Capacitor* capacitor = new Capacitor(name, node1, node2, Coefficient_management(capacitance));
                elements.push_back(capacitor);
                success = true;
            }catch(const exception& e) {
                cout << e.what() << endl;
            }
        }

        void deleteCapacitor(string name) {
            bool found = false;
            for(int i=0; i < elements.size(); i++) {
                if(elements[i]->getName() == name) {
                    elements.erase(elements.begin() + i);
                    found = true;
                }
            }
            if(!found) {
                throw Error("Error: Cannot delete capacitor; component not found");
            }
        }

        void addInductor(string name,string n1,string n2,string inductance) {
            try {
                for(auto e : elements) {
                    if(e->getName() == name) {
                        throw Error("Error: Inductor "+ name +" already exists in the circuit");
                    }
                }
                Node* node1 = new Node(n1,0);
                Node* node2 = new Node(n2,0);
                nodes.push_back(node1 );
                nodes.push_back(node2);
                Inductor* inductor = new Inductor(name, node1, node2, Coefficient_management(inductance));
                elements.push_back(inductor);
                success = true;
            }catch(const exception& e) {
                cout << e.what() << endl;
            }
        }

        void deleteInductor(string name) {
            bool found = false;
            for(int i=0; i < elements.size(); i++) {
                if(elements[i]->getName() == name) {
                    elements.erase(elements.begin() + i);
                    found = true;
                }
            }
            if(!found) {
                throw Error("Error: Cannot delete inductor; component not found");
            }
        }

        void addVoltageSource(string name,string n1,string n2,string voltage) {
            try {
                for(auto e : elements) {
                    if(e->getName() == name) {
                        throw Error("Error: Voltage Source "+ name +" already exists in the circuit");
                    }
                }
                Node* node1 = new Node(n1,0);
                Node* node2 = new Node(n2,0);
                nodes.push_back(node1);
                nodes.push_back(node2);
                VoltageSource* voltageSource = new VoltageSource(name, node1, node2, Coefficient_management(voltage));
                elements.push_back(voltageSource);
                success = true;
            }catch(const exception& e) {
                cout << e.what() << endl;
            }
        }

        void addCurrentSource(string name,string n1,string n2,string current) {
            try {
                for(auto e : elements) {
                    if(e->getName() == name) {
                        throw Error("Error: Current Source "+ name +" already exists in the circuit");
                    }
                }
                Node* node1 = new Node(n1,0);
                Node* node2 = new Node(n2,0);
                nodes.push_back(node1);
                nodes.push_back(node2);
                CurrentSource* currentSource = new CurrentSource(name, node1, node2, Coefficient_management(current));
                elements.push_back(currentSource);
                success = true;
            }catch(const exception& e) {
                cout << e.what() << endl;
            }
        }


        vector<double> gaussianSolve(vector<vector<double>>& A, vector<double>& b) {
            int n = A.size();

            for (int i = 0; i < n; ++i) {
                // Partial pivoting
                int maxRow = i;
                for (int k = i + 1; k < n; ++k) {
                    if (abs(A[k][i]) > abs(A[maxRow][i])) {
                        maxRow = k;
                    }
                }
                swap(A[i], A[maxRow]);
                swap(b[i], b[maxRow]);

                // Eliminate
                for (int k = i + 1; k < n; ++k) {
                    double factor = A[k][i] / A[i][i];
                    for (int j = i; j < n; ++j) {
                        A[k][j] -= factor * A[i][j];
                    }
                    b[k] -= factor * b[i];
                }
            }

            // Back-substitution
            vector<double> x(n);
            for (int i = n - 1; i >= 0; --i) {
                x[i] = b[i];
                for (int j = i + 1; j < n; ++j) {
                    x[i] -= A[i][j] * x[j];
                }
                x[i] /= A[i][i];
            }

            return x;
        }


        vector<vector<double>> LUDecompose(const vector<vector<double>>& A) {
            int n = A.size();
            vector<vector<double>> LU = A;

            for (int i = 0; i < n; ++i) {
                for (int j = i; j < n; ++j) {
                    for (int k = 0; k < i; ++k)
                        LU[i][j] -= LU[i][k] * LU[k][j];
                }

                for (int j = i + 1; j < n; ++j) {
                    for (int k = 0; k < i; ++k)
                        LU[j][i] -= LU[j][k] * LU[k][i];
                    LU[j][i] /= LU[i][i];
                }
            }

            return LU;
        }

        vector<double> LUSolve(const vector<vector<double>>& LU, const vector<double>& b) {
            int n = LU.size();
            vector<double> y(n), x(n);

            // Solve Ly = b (forward substitution)
            for (int i = 0; i < n; ++i) {
                y[i] = b[i];
                for (int k = 0; k < i; ++k)
                    y[i] -= LU[i][k] * y[k];
            }

            // Solve Ux = y (backward substitution)
            for (int i = n - 1; i >= 0; --i) {
                x[i] = y[i];
                for (int k = i + 1; k < n; ++k)
                    x[i] -= LU[i][k] * x[k];
                x[i] /= LU[i][i];
            }

            return x;
        }


        void runTransient(double Tstep, double Tstop, double Tstart = 0.0, double Tmaxstep = -1, const vector<string>& printItems = {}) {
            double time = Tstart + Tstep;
            unordered_map<string, double> capacitorVoltages;
            unordered_map<string, double> inductorCurrents;

            simulationResults.clear();

            while (time <= Tstop + 1e-9) {
                // 1. Create working matrices
                unordered_map<string, int> nodeIndex;
                unordered_map<string, int> voltageSourceIndex;

                int nodeCount = 0;
                for (auto node : nodes) {
                    string name = node->getName();
                    if (name != "0" && name != "GND" && nodeIndex.find(name) == nodeIndex.end()) {
                        nodeIndex[name] = nodeCount++;
                    }
                }

                int n = nodeCount;
                int m = 0;

                for (auto el : elements) {
                    if (el->getType() == "VoltageSource") {
                        voltageSourceIndex[el->getName()] = m++;
                    }
                }

                vector<vector<double>> G(n, vector<double>(n, 0));
                vector<vector<double>> B(n, vector<double>(m, 0));
                vector<vector<double>> C(m, vector<double>(n, 0));
                vector<vector<double>> D(m, vector<double>(m, 0));
                vector<double> J(n, 0);
                vector<double> E(m, 0);

                // 2. Stamp each element
                for (auto el : elements) {
                    string t = el->getType();
                    string n1 = el->getNode1()->getName();
                    string n2 = el->getNode2()->getName();
                    int i = (n1 != "0" && n1 != "GND") ? nodeIndex[n1] : -1;
                    int j = (n2 != "0" && n2 != "GND") ? nodeIndex[n2] : -1;

                    if (t == "Resistor") {
                        double R = el->getValue();
                        double y = 1.0 / R;
                        if (i != -1) G[i][i] += y;
                        if (j != -1) G[j][j] += y;
                        if (i != -1 && j != -1) {
                            G[i][j] -= y;
                            G[j][i] -= y;
                        }
                    } else if (t == "Capacitor") {
                        double Cval = el->getValue();
                        double g = Cval / Tstep;
                        double v_prev = capacitorVoltages[el->getName()];
                        double Ieq = g * v_prev;

                        if (i != -1) G[i][i] += g;
                        if (j != -1) G[j][j] += g;
                        if (i != -1 && j != -1) {
                            G[i][j] -= g;
                            G[j][i] -= g;
                        }
                        if (i != -1) J[i] += Ieq;
                        if (j != -1) J[j] -= Ieq;
                    } else if (t == "Inductor") {
                        double Lval = el->getValue();
                        double g = Tstep / Lval;
                        double i_prev = inductorCurrents[el->getName()];
                        double Veq = i_prev * g;

                        if (i != -1) G[i][i] += g;
                        if (j != -1) G[j][j] += g;
                        if (i != -1 && j != -1) {
                            G[i][j] -= g;
                            G[j][i] -= g;
                        }
                        if (i != -1) J[i] += Veq;
                        if (j != -1) J[j] -= Veq;
                    } else if (t == "CurrentSource") {
                        double Ival = el->getValue();
                        if (i != -1) J[i] -= Ival;
                        if (j != -1) J[j] += Ival;
                    } else if (t == "VoltageSource") {
                        int k = voltageSourceIndex[el->getName()];
                        double V = el->getValue();
                        if (i != -1) {
                            B[i][k] = 1;
                            C[k][i] = 1;
                        }
                        if (j != -1) {
                            B[j][k] = -1;
                            C[k][j] = -1;
                        }
                        E[k] = V;
                    }
                }

                // 3. Form A and Z
                int size = n + m;
                vector<vector<double>> A(size, vector<double>(size, 0));
                vector<double> Z(size, 0);

                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j)
                        A[i][j] = G[i][j];

                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < m; ++j)
                        A[i][n + j] = B[i][j];

                for (int i = 0; i < m; ++i)
                    for (int j = 0; j < n; ++j)
                        A[n + i][j] = C[i][j];

                for (int i = 0; i < m; ++i)
                    for (int j = 0; j < m; ++j)
                        A[n + i][n + j] = D[i][j];

                for (int i = 0; i < n; ++i)
                    Z[i] = J[i];
                for (int i = 0; i < m; ++i)
                    Z[n + i] = E[i];

                // 4. Solve
                vector<vector<double>> LU = LUDecompose(A);
                vector<double> X = LUSolve(LU, Z);

                // 5. Store results
                TimePointResult result;
                result.time = time;
                for (auto it = nodeIndex.begin(); it != nodeIndex.end(); ++it) {
                    string name = it->first;
                    int idx = it->second;
                    result.nodeVoltages[name] = X[idx];
                }

                for (const auto& el : elements) {
                    if (el->getType() == "VoltageSource") {
                        string name = el->getName();
                        int idx = voltageSourceIndex[name];
                        result.elementCurrents[name] = X[n + idx];
                    }
                }

                for (const auto& el : elements) {
                    string name = el->getName();
                    string type = el->getType();
                    string n1 = el->getNode1()->getName();
                    string n2 = el->getNode2()->getName();

                    int i = (n1 != "0" && n1 != "GND") ? nodeIndex[n1] : -1;
                    int j = (n2 != "0" && n2 != "GND") ? nodeIndex[n2] : -1;

                    double v1 = (i != -1) ? X[i] : 0.0;
                    double v2 = (j != -1) ? X[j] : 0.0;
                    double voltageAcross = v1 - v2;

                    if (type == "Resistor") {
                        double R = el->getValue();
                        double current = voltageAcross / R;
                        result.elementCurrents[name] = current;

                    } else if (type == "Capacitor") {
                        double C = el->getValue();
                        double v_prev = capacitorVoltages[name];
                        double current = C * (voltageAcross - v_prev) / Tstep;
                        result.elementCurrents[name] = current;

                    } else if (type == "Inductor") {
                        result.elementCurrents[name] = inductorCurrents[name];
                    }
                }

                simulationResults.push_back(result);

                // 6. Update history for next time step
                for (auto el : elements) {
                    if (el->getType() == "Capacitor") {
                        string name = el->getName();
                        int i = (el->getNode1()->getName() != "0" && el->getNode1()->getName() != "GND") ? nodeIndex[el->getNode1()->getName()] : -1;
                        int j = (el->getNode2()->getName() != "0" && el->getNode2()->getName() != "GND") ? nodeIndex[el->getNode2()->getName()] : -1;
                        double v1 = (i != -1) ? X[i] : 0;
                        double v2 = (j != -1) ? X[j] : 0;
                        capacitorVoltages[name] = v1 - v2;
                    } else if (el->getType() == "Inductor") {
                        string name = el->getName();
                        int i = (el->getNode1()->getName() != "0" && el->getNode1()->getName() != "GND") ? nodeIndex[el->getNode1()->getName()] : -1;
                        int j = (el->getNode2()->getName() != "0" && el->getNode2()->getName() != "GND") ? nodeIndex[el->getNode2()->getName()] : -1;
                        double v1 = (i != -1) ? X[i] : 0;
                        double v2 = (j != -1) ? X[j] : 0;
                        inductorCurrents[name] += (v1 - v2) * (Tstep / el->getValue());
                    }
                }

                time += Tstep;
            }


            // 7. Print simulation output (filtered)
            vector<string> voltageNodes;
            vector<string> currentElements;

            for (const auto& item : printItems) {
                if (item.size() > 2 && item.substr(0, 2) == "V(" && item.back() == ')') {
                    string nodeName = item.substr(2, item.size() - 3);
                    voltageNodes.push_back(nodeName);
                } else if (item.size() > 2 && item.substr(0, 2) == "I(" && item.back() == ')') {
                    string elementName = item.substr(2, item.size() - 3);
                    currentElements.push_back(elementName);
                }
            }

            for (const auto& n : voltageNodes) cout << "V(" << n << ")\t";
            for (const auto& e : currentElements) cout << "I(" << e << ")\t";
            cout << endl;

            for (const auto& res : simulationResults) {
                for (const auto& n : voltageNodes) {
                    if (res.nodeVoltages.count(n))
                        cout << res.nodeVoltages.at(n) << "\t";
                    else
                        cout << "N/A\t";
                }
                for (const auto& e : currentElements) {
                    if (res.elementCurrents.count(e))
                        cout << res.elementCurrents.at(e) << "\t";
                    else
                        cout << "N/A\t";
                }
                cout << endl;
            }
        }

        void clear() {
            elements.clear();
            nodes.clear();
        }
    };
}

using namespace Controller;
//---------------------------------------------------------------------------------------------------------------------------------------
namespace fs = filesystem;

class FileHandler {
private:
    string folder = "schematics";

public:
    FileHandler() {
        if (!fs::exists(folder))
            fs::create_directory(folder);
    }

    vector<string> listSchematics() {
        vector<string> files;
        int idx = 1;
        for (const auto &entry : fs::directory_iterator(folder)) {
            if (entry.is_regular_file()) {
                string name = entry.path().filename().string();
                files.push_back(name);
                cout << idx++ << "- " << name << endl;
            }
        }
        return files;
    }

    bool isValidIndex(int index, int max) {
        return index >= 1 && index <= max;
    }

    bool loadSchematic(const string& filename, Circuit* circuit, vector<string> &netlistLines) {
        ifstream infile(folder + "/" + filename);
        if (!infile.is_open()) {
            cout << "-Error: Cannot open schematic file." << endl;
            return false;
        }

        string line;
        while (getline(infile, line)) {
            if (line.empty()) continue;

            cout << line << endl; // Print if needed
            netlistLines.push_back(line);

            istringstream iss(line);
            string name, n1, n2, value;
            iss >> name >> n1 >> n2 >> value;

            try {
                if (name[0] == 'R')
                    circuit->addResistor(name, n1, n2, value);
                else if (name[0] == 'C')
                    circuit->addCapacitor(name, n1, n2, value);
                else if (name[0] == 'L')
                    circuit->addInductor(name, n1, n2, value);
                else if (name[0] == 'V')
                    circuit->addVoltageSource(name, n1, n2, value);
                else if (name[0] == 'I')
                    circuit->addCurrentSource(name, n1, n2, value);
                else
                    throw runtime_error("Unknown element type: " + name);
            } catch (const exception& e) {
                cout << "-Error: " << e.what() << endl;
            }
        }

        return true;
    }

    void saveSchematic(const string& filename, const vector<string>& netlistLines) {
        ofstream outfile(folder + "/" + filename);
        for (const auto& line : netlistLines) {
            outfile << line << endl;
        }
    }
};


//---------------------------------------------------------------------------------------------------------------------------------------
namespace View {
    class view {
    private:
        Circuit *circuit;

        enum class MenuState {
            MAIN_MENU,
            SCHEMATIC_SELECTED
        };
        MenuState menu = MenuState::MAIN_MENU;
    public:
        view(Circuit *c) : circuit(c) {}

        void run() {
            FileHandler fileHandler;
            string activeFile;
            vector<string> netlistLines;

            while (true) {
                success = false;

                try {
                    if (menu == MenuState::MAIN_MENU) {
                        cout << "\n- Main Menu:\n";
                        cout << "1. show existing schematics\n";
                        cout << "2. NewFile <file_path>\n";
                        cout << "Type 'quit' to exit\n>>> ";

                        string command;
                        getline(cin, command);

                        if (command == "quit")
                            break;

                        if (command == "1") {
                            auto files = fileHandler.listSchematics();
                            if (files.empty()) {
                                cout << "No schematics found.\n";
                                continue;
                            }

                            cout << "- choose existing schematic:\n>>> ";
                            string selection;
                            getline(cin, selection);

                            try {
                                int index = stoi(selection);
                                if (!fileHandler.isValidIndex(index, files.size()))
                                    throw invalid_argument("Invalid index");

                                activeFile = files[index - 1];
                                circuit->clear();
                                netlistLines.clear();
                                if (fileHandler.loadSchematic(activeFile, circuit, netlistLines)) {
                                    menu = MenuState::SCHEMATIC_SELECTED;
                                    continue;
                                }
                            } catch (...) {
                                cout << "-Error: Inappropriate input\n";
                            }

                        } else if (command.find("NewFile ") == 0) {
                            string filePath = command.substr(8);

                            ofstream newFile(filePath); // This will create the file if it doesn't exist
                            if (!newFile.is_open()) {
                                cout << "-Error: Cannot create file\n";
                                continue;
                            }


                            cout << "new file created successfully in " << filePath << endl;
                        } else {
                            cout << "-Error: Inappropriate input\n";
                        }

                    } else if (menu == MenuState::SCHEMATIC_SELECTED) {
                        cout << "\nEnter schematic command or 'return' to go back to menu:\n>>> ";
                        string input;
                        getline(cin, input);

                        if (input == "return") {
                            menu = MenuState::MAIN_MENU;
                            continue;
                        }

                        regex transient(".TRAN (\\S+) (\\S+) (\\S+) (\\S+)");
                        smatch match1;

                        if (regex_match(input, match1, transient)) {
                            double Tstep = stod(match1[1]);
                            double Tstop = stod(match1[2]);
                            string voltage = match1[3];
                            string current = match1[4];

                            vector<string> printItems = { voltage, current };
                            circuit->runTransient(Tstep, Tstop, 0, -1, printItems);

                            bool analysis = false;
                            for (auto l : netlistLines) {
                                if (l.substr(0, 5) == ".TRAN"){
                                    analysis = true;
                                    break;
                                }
                            }

                            if (!analysis)
                                netlistLines.push_back(".TRAN " + match1[1].str() + " " + match1[2].str());
                            fileHandler.saveSchematic(activeFile, netlistLines);
                            continue;
                        }

                        regex add("add (\\S+) (\\S+) (\\S+) (\\S+)");
                        regex del("delete (\\S+)");
                        smatch match;

                        if (regex_search(input, match, add)) {
                            try {
                                string id = match[1];
                                string n1 = match[2], n2 = match[3], val = match[4];

                                string fileInput = id + " " + n1 + " " + n2 + " " + val;

                                if (id[0] == 'R'){
                                    circuit->addResistor(id, n1, n2, val);
                                } else if (id[0] == 'C') {
                                    circuit->addCapacitor(id, n1, n2, val);
                                } else if (id[0] == 'L'){
                                    circuit->addInductor(id, n1, n2, val);
                                } else if (id[0] == 'V'){
                                    circuit->addVoltageSource(id, n1, n2, val);
                                } else if (id[0] == 'I'){
                                    circuit->addCurrentSource(id, n1, n2, val);
                                } else
                                    throw Error("Error: Element type not found.");


                                if (success)
                                    netlistLines.push_back(fileInput);
                                fileHandler.saveSchematic(activeFile, netlistLines);


                            } catch (const exception& e) {
                                cout << e.what() << endl;
                            }

                        } else if (regex_search(input, match, del)) {
                            try {
                                string id = match[1];
                                if (id[0] == 'R') {
                                    circuit->deleteResistor(id);
                                    netlistLines.erase(remove_if(netlistLines.begin(), netlistLines.end(),
                                                                 [&](const string& line) { return line.find(id) != string::npos; }),
                                                       netlistLines.end());
                                } else if (id[0] == 'C') {
                                    circuit->deleteCapacitor(id);
                                    netlistLines.erase(remove_if(netlistLines.begin(), netlistLines.end(),
                                                                 [&](const string& line) { return line.find(id) != string::npos; }),
                                                       netlistLines.end());
                                } else if (id[0] == 'L') {
                                    circuit->deleteInductor(id);
                                    netlistLines.erase(remove_if(netlistLines.begin(), netlistLines.end(),
                                                                 [&](const string& line) { return line.find(id) != string::npos; }),
                                                       netlistLines.end());
                                } else
                                    throw Error("Error: Element not found.");


                                fileHandler.saveSchematic(activeFile, netlistLines);
                            } catch (const exception& e) {
                                cout << e.what() << endl;
                            }

                        } else {
                            cout << "-Error: Syntax error\n";
                        }
                    }

                } catch (const exception& e) {
                    cout << e.what() << endl;
                }
            }

        }
    };
}


using namespace View;

int main() {
    Circuit circuits;
    view Views(&circuits);
    Views.run();
    return 0;
}



//add R1 GND N2 2
//add R2 GND N2 2
//add V2 GND N2 4
//.TRAN 0.01 5 V(N2) I(R1)