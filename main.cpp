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
        string name = "";
        double voltage = 0;
        bool isGround = false;


    public:
        Node(string name, double voltage) : name(name), voltage(voltage) {}

        string getName() { return name; }

        double getVoltage() {
            return voltage;
        }

        void setVoltage(double v) {
            voltage = v;
        }

        bool getIsGround (){
            return isGround;
        }

        void setIsGround (bool state){
            isGround=state;
        }

    };


    //---------------------------------------------------------------------------------------------------------------------------------------
    class Element {
    private:
        Node *node1;
        Node *node2;
        string name = "";
        string type = "";
        double voltage = 0;
        double current = 0;

    public:
        Element(string name, string type, Node *node1, Node *node2)
                : name(name), type(type), node1(node1), node2(node2) {}

        string getName() const {
            return name;
        }

        string getType() const {
            return type;
        }

        Node *getNode1() {
            return node1;
        }

        Node *getNode2() {
            return node2;
        }

        double getVoltage(){
            return node1->getVoltage() - node2->getVoltage();
        }

        virtual double getValue() const {}

        virtual void setValue(double r) {}

        void setVoltage(double voltage){ this->voltage= voltage;}

        void setCurrent(double current){ this->current= current;}

        double getCurrent(){return current;}

        void setType(const string& t) { type = t; }

        virtual ~Element() = default;
    };


    class Resistor : public Element {
    private:
        double resistance = 0;

    public:
        Resistor(string name, Node *node1, Node *node2, double resistance)
                : Element(name, "Resistor", node1, node2), resistance(resistance) {}

        double getValue() const override {
            return resistance;
        }
    };


    class Inductor : public Element {
    private:
        double inductance = 0;
        double previousCurrent = 0.0;

    public:
        Inductor(string name, Node *node1, Node *node2, double inductance)
                : Element(name, "Inductor", node1, node2), inductance(inductance) {}

        double getValue() const override {
            return inductance;
        }

        double getPreviousCurrent (){
            return previousCurrent;
        }
        void setPreviousCurrent (double d){
            previousCurrent = d;
        }
    };


    class Capacitor : public Element {
    private:
        double capacitance = 0;
        double previousVoltage = 0.0;

    public:
        Capacitor(string name, Node *node1, Node *node2, double capacitance)
                : Element(name, "Capacitor", node1, node2), capacitance(capacitance) {}

        double getValue() const override {
            return capacitance;
        }

        double getPreviousVoltage(){
            return previousVoltage;
        }
        void setPreviousVoltage (double d){
            previousVoltage = d;
        }
    };


    class CurrentSource : public Element {
    private:
        double current = 1;

    public:
        CurrentSource(string name, Node *node1, Node *node2, double current)
                : Element(name, "CurrentSource", node1, node2), current(current) {}

        double getValue() const override { return current; }
        void setValue(double c) override { current = c; }
        virtual void setValueAtTime(double t){
        }
    };


    class VoltageSource : public Element {
    private:
        double voltage = 1;

    public:
        VoltageSource(string name, Node *node1, Node *node2, double voltage)
                : Element(name, "VoltageSource", node1, node2), voltage(voltage) {}

        double getValue() const override { return voltage; }
        void setValue(double c) override { voltage = c; }
        virtual void setValueAtTime(double t){
        }
    };

    class Diode : public Element {
    private:
        string model;
        double vOn;
        bool isOn;
        bool assumedOn = false;
    public:
        double getValue() const override { return vOn; }

        void setValue(double r) override {}

        Diode(string name = "", string mode = "", Node *n = nullptr, Node *p = nullptr, double v = 0)
                : Element("Diode" + mode, name, n, p), vOn(v), isOn(false) {}

        void updateState() {
            double v = getVoltage();
            isOn = (v >= vOn);
        }

        bool getIsOn() const { return isOn; }

        void assumeState(bool on) {
            this->assumedOn = on;
        }

        bool isAssumedOn() const { return assumedOn; }
    };

    class SinVoltageSource : public VoltageSource {
    private:
        double Voffset;
        double Vamplitude;
        double Frequency;
    public:
        SinVoltageSource(
                string name = "",
                Node *n = nullptr,
                Node *p = nullptr,
                double offset = 0.0,
                double amplitude = 0.0,
                double freq = 0.0
        )
                : VoltageSource(name, n, p, 0.0),
                  Voffset(offset), Vamplitude(amplitude), Frequency(freq) {}

        double getOffset() const { return Voffset; }

        double getAmplitude() const { return Vamplitude; }

        double getFrequency() const { return Frequency; }

        void setValueAtTime(double t) {
            setValue(Voffset + Vamplitude * sin(2 * M_PI * Frequency * t));
        }

        void setOffset(double offset) { Voffset = offset; }

        void setAmplitude(double amplitude) { Vamplitude = amplitude; }

        void setFrequency(double freq) { Frequency = freq; }
    };

    class SinCurrentSource : public CurrentSource {
    private:
        double Voffset;
        double Vamplitude;
        double Frequency;
    public:
        SinCurrentSource(
                string name = "",
                Node *n = nullptr,
                Node *p = nullptr,
                double offset = 0.0,
                double amplitude = 0.0,
                double freq = 0.0
        )
                : CurrentSource(name, n, p, 0.0),
                  Voffset(offset), Vamplitude(amplitude), Frequency(freq) {}

        double getOffset() const { return Voffset; }

        double getAmplitude() const { return Vamplitude; }

        double getFrequency() const { return Frequency; }

        void setValueAtTime(double t) {
            setValue(Voffset + Vamplitude * sin(2 * M_PI * Frequency * t));
        }

        void setOffset(double offset) { Voffset = offset; }

        void setAmplitude(double amplitude) { Vamplitude = amplitude; }

        void setFrequency(double freq) { Frequency = freq; }
    };

    class PulseVoltageSource : public VoltageSource {
    private:
        double Vinitial;
        double Von;
        double Tdelay;
        double Trise;
        double Tfall;
        double Ton;
        double Tperiod;
        double Ncycles;
    public:
        PulseVoltageSource(
                const string &name = "",
                Node *n = nullptr,
                Node *p = nullptr,
                double Vinitial = 0.0,
                double Von = 0.0,
                double Tdelay = 0.0,
                double Trise = 0.0,
                double Tfall = 0.0,
                double Ton = 0.0,
                double Tperiod = 0.0,
                double Ncycles = 0.0
        )
                : VoltageSource(name, n, p, 0.0),
                  Vinitial(Vinitial), Von(Von), Tdelay(Tdelay),
                  Trise(Trise), Tfall(Tfall), Ton(Ton),
                  Tperiod(Tperiod), Ncycles(Ncycles) {}

        double getVinitial() const { return Vinitial; }

        double getVon() const { return Von; }

        double getTdelay() const { return Tdelay; }

        double getTrise() const { return Trise; }

        double getTfall() const { return Tfall; }

        double getTon() const { return Ton; }

        double getTperiod() const { return Tperiod; }

        double getNcycles() const { return Ncycles; }

        void setValueAtTime(double t) {
            if (t < Tdelay) setValue(Vinitial);
            double localTime = t - Tdelay;
            double cycleTime = fmod(localTime, Tperiod);
            if (Ncycles != 0 && localTime > Ncycles * Tperiod) {
                setValue(Vinitial);
            }
            if (cycleTime < Trise) {
                setValue(Vinitial + (Von - Vinitial) * (cycleTime / Trise));
            } else if (cycleTime < Trise + Ton) {
                setValue(Von);
            } else if (cycleTime < Trise + Ton + Tfall) {
                setValue(Von - (Von - Vinitial) * ((cycleTime - Trise - Ton) / Tfall));
            } else {
                setValue(Vinitial);
            }
        }

        void setVinitial(double v) { Vinitial = v; }

        void setVon(double v) { Von = v; }

        void setTdelay(double t) { Tdelay = t; }

        void setTrise(double t) { Trise = t; }

        void setTfall(double t) { Tfall = t; }

        void setTon(double t) { Ton = t; }

        void setTperiod(double t) { Tperiod = t; }

        void setNcycles(double n) { Ncycles = n; }
    };

    class PulseCurrentSource : public CurrentSource {
    private:
        double Iinitial;
        double Ion;
        double Tdelay;
        double Trise;
        double Tfall;
        double Ton;
        double Tperiod;
        double Ncycles;
    public:
        PulseCurrentSource(
                const string &name = "",
                Node *n = nullptr,
                Node *p = nullptr,
                double Iinitial = 0.0,
                double Ion = 0.0,
                double Tdelay = 0.0,
                double Trise = 0.0,
                double Tfall = 0.0,
                double Ton = 0.0,
                double Tperiod = 0.0,
                double Ncycles = 0.0
        )
                : CurrentSource(name, n, p, 0.0),
                  Iinitial(Iinitial), Ion(Ion), Tdelay(Tdelay),
                  Trise(Trise), Tfall(Tfall), Ton(Ton),
                  Tperiod(Tperiod), Ncycles(Ncycles) {}

        double getIinitial() const { return Iinitial; }

        double getIon() const { return Ion; }

        double getTdelay() const { return Tdelay; }

        double getTrise() const { return Trise; }

        double getTfall() const { return Tfall; }

        double getTon() const { return Ton; }

        double getTperiod() const { return Tperiod; }

        double getNcycles() const { return Ncycles; }

        void setValueAtTime(double t) {
            if (t < Tdelay) setValue(Iinitial);
            double localTime = t - Tdelay;
            double cycleTime = fmod(localTime, Tperiod);
            if (Ncycles != 0 && localTime > Ncycles * Tperiod) {
                setValue(Iinitial);
            }
            if (cycleTime < Trise) {
                setValue(Iinitial + (Ion - Iinitial) * (cycleTime / Trise));
            } else if (cycleTime < Trise + Ton) {
                setValue(Ion);
            } else if (cycleTime < Trise + Ton + Tfall) {
                setValue(Ion - (Ion - Iinitial) * ((cycleTime - Trise - Ton) / Tfall));
            } else {
                setValue(Iinitial);
            }
        }

        void setIinitial(double i) { Iinitial = i; }

        void setIon(double i) { Ion = i; }

        void setTdelay(double t) { Tdelay = t; }

        void setTrise(double t) { Trise = t; }

        void setTfall(double t) { Tfall = t; }

        void setTon(double t) { Ton = t; }

        void setTperiod(double t) { Tperiod = t; }

        void setNcycles(double n) { Ncycles = n; }
    };

    class deltaVoltageSource : public VoltageSource {
    private:
        double tPulse;
        double area;
        double step = 0.001;
    public:
        deltaVoltageSource(const std::string &name, Node *p, Node *n,
                           double tPulse, double step, double area = 1.0)
                : VoltageSource(name, n, p, 0.0), tPulse(tPulse), area(area) {}

        void setValueAtTime(double t) override {
            if (std::fabs(t - tPulse) < 1e-12) {
                setVoltage(area / step);
            } else {
                setVoltage(0.0);
            }
        }

        void setStep(double ste) {
            step = ste;
        }
    };

    class deltaCurrentSource : public CurrentSource {
    private:
        double tPulse;
        double area;
        double step;
    public:
        deltaCurrentSource(const std::string &name, Node *p, Node *n,
                           double tPulse, double step, double area = 1.0)
                : CurrentSource(name, n, p, 0.0), tPulse(tPulse), area(area) {}

        void setValueAtTime(double t) override {
            if (std::fabs(t - tPulse) < 1e-12) {
                setCurrent(area / step);
            } else {
                setCurrent(0.0);
            }
        }

        void setStep(double ste) {
            step = ste;
        }
    };

    class CCCS : public CurrentSource {
    private:
        double gain;
        Element *element;
    public:
        CCCS(std::string name = "", Node *n = nullptr, Node *p = nullptr, double gain = 1.0, Element *e = nullptr)
                : CurrentSource(name, n, p, 0.0), gain(gain), element(e) {
            setType("CCCS");
        }

        double getGain() {
            return gain;
        }

        Element *getControlSourceName() {
            return element;
        }

        void setValueAtTime(double t = 0) {
            CurrentSource::setValue(element->getCurrent());
        }
    };

    class CCVS : public VoltageSource {
    private:
        double gain;
        Element *element;

    public:
        CCVS(std::string name = "", Node *n = nullptr, Node *p = nullptr, double gain = 1.0, Element *e = nullptr)
                : VoltageSource(name, n, p, 0.0), gain(gain), element(e) {
            setType("CCVS");
        }

        double getGain() const {
            return gain;
        }

        Element *getControlSourceName() const {
            return element;
        }

        void setValueAtTime(double t = 0) {
            if (element != nullptr) {
                setValue(gain * element->getCurrent());
            }
        }
    };

    class VCCS : public CurrentSource {
    private:
        double gain;
        Node *ctrlPositive;
        Node *ctrlNegative;
    public:
        VCCS(std::string name = "", Node *n = nullptr, Node *p = nullptr,
             double gain = 1.0, Node *ctrlPos = nullptr, Node *ctrlNeg = nullptr)
                : CurrentSource(name, n, p, 0.0), gain(gain), ctrlPositive(ctrlPos), ctrlNegative(ctrlNeg) {
            setType("VCCS");
        }

        double getGain() const {
            return gain;
        }

        void setValueAtTime(double t = 0) {
            if (ctrlPositive && ctrlNegative) {
                double vctrl = ctrlPositive->getVoltage() - ctrlNegative->getVoltage();
                setValue(gain * vctrl);
            }
        }

        Node *getControlNodeP() const { return ctrlPositive; }

        Node *getControlNodeN() const { return ctrlNegative; }
    };

    class VCVS : public VoltageSource {
    private:
        double gain;
        Node *ctrlPositive;
        Node *ctrlNegative;

    public:
        VCVS(std::string name = "", Node *n = nullptr, Node *p = nullptr,
             double gain = 1.0, Node *ctrlPos = nullptr, Node *ctrlNeg = nullptr)
                : VoltageSource(name, n, p, 0.0), gain(gain), ctrlPositive(ctrlPos), ctrlNegative(ctrlNeg) {
            setType("VCVS");
        }

        double getGain() const {
            return gain;
        }

        void setValueAtTime(double t = 0) {
            if (ctrlPositive && ctrlNegative) {
                double vctrl = ctrlPositive->getVoltage() - ctrlNegative->getVoltage();
                setValue(gain * vctrl);
            }
        }

        Node *getControlNodeP() const { return ctrlPositive; }

        Node *getControlNodeN() const { return ctrlNegative; }
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

        double unitHandler(string input, string type) {
            regex pattern(R"(^\s*([+-]?[\d]*\.?[\d]+(?:[eE][+-]?[\d]+)?)\s*(Meg|[numkMG])?\s*$)");
            smatch match;
            if (regex_match(input, match, pattern)) {
                double value = stod(match[1]);
                string unit = match[2].str();
                double multiplier = 1.0;
                if (unit == "n") multiplier = 1e-9;
                else if (unit == "u") multiplier = 1e-6;
                else if (unit == "m") multiplier = 1e-3;
                else if (unit == "k") multiplier = 1e3;
                else if (unit == "M" || unit == "Meg") multiplier = 1e6;
                else if (unit == "G") multiplier = 1e9;
                if (value <= 0) {
                    throw Error("Error: Negative or zero value for a component is invalid.");
                }
                return value * multiplier;
            } else {
                throw Error("Error: Negative or zero value for a component is invalid.");
            }
        }

        double unitHandler2(string input, string type="time") {
            regex pattern(R"(^\s*([+-]?[\d]*\.?[\d]+(?:[eE][+-]?[\d]+)?)\s*(Meg|[numkMG])?\s*$)");
            smatch match;
            if (regex_match(input, match, pattern)) {
                double value = stod(match[1]);
                string unit = match[2].str();
                double multiplier = 1.0;
                if (unit == "n") multiplier = 1e-9;
                else if (unit == "u") multiplier = 1e-6;
                else if (unit == "m") multiplier = 1e-3;
                else if (unit == "k") multiplier = 1e3;
                else if (unit == "M" || unit == "Meg") multiplier = 1e6;
                else if (unit == "G") multiplier = 1e9;
                if (value < 0) {
                    throw Error("Error: Negative or zero value for a component is invalid.");
                }
                return value * multiplier;
            } else {
                throw Error("Error: Negative or zero value for a component is invalid.");
            }
        }

        double unitHandler3(string input, string type) {
            regex pattern(R"(^\s*([+-]?[\d]*\.?[\d]+(?:[eE][+-]?[\d]+)?)\s*(Meg|[numkMG])?\s*$)");
            smatch match;
            if (regex_match(input, match, pattern)) {
                double value = stod(match[1]);
                string unit = match[2].str();
                double multiplier = 1.0;
                if (unit == "n") multiplier = 1e-9;
                else if (unit == "u") multiplier = 1e-6;
                else if (unit == "m") multiplier = 1e-3;
                else if (unit == "k") multiplier = 1e3;
                else if (unit == "M" || unit == "Meg") multiplier = 1e6;
                else if (unit == "G") multiplier = 1e9;
                return value * multiplier;
            } else {
                throw Error("Error: Negative or zero value for a component is invalid.");
            }
        }

        Node* findNode(string name){
            for (Node* &i:nodes) {
                if(i->getName()==name){return i;}
            }
            return nullptr;
        }

        Element* findElement (string name){
            for (int i = 0; i < elements.size(); ++i) {
                if(elements[i]->getName()==name) {
                    return elements[i];
                }
            }
            return nullptr;
        }

        void addResistor(string name, string n1, string n2, string resistance) {
            try {
                for (auto e : elements) {
                    if (e->getName() == name) {
                        throw Error("Error: Resistor " + name + " already exists in the circuit");
                    }
                }

                // Find or create node1
                Node* node1 = new Node(n1, 0);
                auto it1 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n1; });
                if (it1 == nodes.end()) {
                    nodes.push_back(node1);
                } else {
                    delete node1;
                    node1 = *it1;
                }

                // Find or create node2
                Node* node2 = new Node(n2, 0);
                auto it2 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n2; });
                if (it2 == nodes.end()) {
                    nodes.push_back(node2);
                } else {
                    delete node2;
                    node2 = *it2;
                }

                // Add voltage source element
                Resistor* resistor = new Resistor(name, node1, node2, Coefficient_management(resistance));
                elements.push_back(resistor);

                success = true;
            } catch (const exception& e) {
                cout << e.what() << endl;
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

        void addCapacitor(string name, string n1, string n2, string capacitance) {
            try {
                for (auto e : elements) {
                    if (e->getName() == name) {
                        throw Error("Error: Capacitor " + name + " already exists in the circuit");
                    }
                }

                // Find or create node1
                Node* node1 = new Node(n1, 0);
                auto it1 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n1; });
                if (it1 == nodes.end()) {
                    nodes.push_back(node1);
                } else {
                    delete node1;
                    node1 = *it1;
                }

                // Find or create node2
                Node* node2 = new Node(n2, 0);
                auto it2 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n2; });
                if (it2 == nodes.end()) {
                    nodes.push_back(node2);
                } else {
                    delete node2;
                    node2 = *it2;
                }

                // Add voltage source element
                Capacitor* capacitor = new Capacitor(name, node1, node2, Coefficient_management(capacitance));
                elements.push_back(capacitor);

                success = true;
            } catch (const exception& e) {
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

        void addInductor(string name, string n1, string n2, string inductance) {
            try {
                for (auto e : elements) {
                    if (e->getName() == name) {
                        throw Error("Error: Inductance " + name + " already exists in the circuit");
                    }
                }

                // Find or create node1
                Node* node1 = new Node(n1, 0);
                auto it1 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n1; });
                if (it1 == nodes.end()) {
                    nodes.push_back(node1);
                } else {
                    delete node1;
                    node1 = *it1;
                }

                // Find or create node2
                Node* node2 = new Node(n2, 0);
                auto it2 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n2; });
                if (it2 == nodes.end()) {
                    nodes.push_back(node2);
                } else {
                    delete node2;
                    node2 = *it2;
                }

                // Add voltage source element
                Inductor* inductor = new Inductor(name, node1, node2, Coefficient_management(inductance));
                elements.push_back(inductor);

                success = true;
            } catch (const exception& e) {
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

        void addVoltageSource(string name, string n1, string n2, string voltage) {
            try {
                for (auto e : elements) {
                    if (e->getName() == name) {
                        throw Error("Error: Voltage Source " + name + " already exists in the circuit");
                    }
                }

                // Find or create node1
                Node* node1 = new Node(n1, 0);
                auto it1 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n1; });
                if (it1 == nodes.end()) {
                    nodes.push_back(node1);
                } else {
                    delete node1;
                    node1 = *it1;
                }

                // Find or create node2
                Node* node2 = new Node(n2, 0);
                auto it2 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n2; });
                if (it2 == nodes.end()) {
                    nodes.push_back(node2);
                } else {
                    delete node2;
                    node2 = *it2;
                }

                // Add voltage source element
                VoltageSource* voltageSource = new VoltageSource(name, node1, node2, Coefficient_management(voltage));
                elements.push_back(voltageSource);

                success = true;
            } catch (const exception& e) {
                cout << e.what() << endl;
            }
        }

        void addCurrentSource(string name, string n1, string n2, string current) {
            try {
                // Check for duplicate element
                for (auto e : elements) {
                    if (e->getName() == name) {
                        throw Error("Error: Current Source " + name + " already exists in the circuit");
                    }
                }

                // Find or create node1
                Node* node1 = new Node(n1, 0);
                auto it1 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n1; });
                if (it1 == nodes.end()) {
                    nodes.push_back(node1);
                } else {
                    delete node1;
                    node1 = *it1;
                }

                // Find or create node2
                Node* node2 = new Node(n2, 0);
                auto it2 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n2; });
                if (it2 == nodes.end()) {
                    nodes.push_back(node2);
                } else {
                    delete node2;
                    node2 = *it2;
                }

                // Add current source element
                CurrentSource* currentSource = new CurrentSource(name, node1, node2, Coefficient_management(current));
                elements.push_back(currentSource);

                success = true;
            } catch (const exception& e) {
                cout << e.what() << endl;
            }
        }


        void addVCVS(string name, string n1, string n2, string cn1, string cn2, string gain) {
            try {
                // Check for duplicate element name
                for (auto e : elements) {
                    if (e->getName() == name) {
                        throw Error("Error: Element " + name + " already exists in the circuit");
                    }
                }

                // Find or create node1
                Node* node1 = new Node(n1, 0);
                auto it1 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n1; });
                if (it1 == nodes.end()) {
                    nodes.push_back(node1);
                } else {
                    delete node1; // Avoid memory leak
                    node1 = *it1;
                }

                // Find or create node2
                Node* node2 = new Node(n2, 0);
                auto it2 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n2; });
                if (it2 == nodes.end()) {
                    nodes.push_back(node2);
                } else {
                    delete node2;
                    node2 = *it2;
                }

                // Control nodes must exist
                auto ctrlIt1 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == cn1; });
                if (ctrlIt1 == nodes.end()) {
                    throw Error("Error: Control node " + cn1 + " not found");
                }

                auto ctrlIt2 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == cn2; });
                if (ctrlIt2 == nodes.end()) {
                    throw Error("Error: Control node " + cn2 + " not found");
                }

                // Create and add VCVS element
                double g = Coefficient_management(gain);
                VCVS* vcvs = new VCVS(name, node2, node1, g, *ctrlIt1, *ctrlIt2);
                elements.push_back(vcvs);

                success = true;
            } catch (const exception& e) {
                cout << e.what() << endl;
            }
        }

        void addVCCS(string name, string n1, string n2, string cn1, string cn2, string gain) {
            try {
                // Check for duplicate element
                for (auto e : elements) {
                    if (e->getName() == name) {
                        throw Error("Error: Element " + name + " already exists in the circuit");
                    }
                }

                // Find or create node1
                Node* node1 = new Node(n1, 0);
                auto it1 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n1; });
                if (it1 == nodes.end()) {
                    nodes.push_back(node1);
                } else {
                    delete node1;
                    node1 = *it1;
                }

                // Find or create node2
                Node* node2 = new Node(n2, 0);
                auto it2 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n2; });
                if (it2 == nodes.end()) {
                    nodes.push_back(node2);
                } else {
                    delete node2;
                    node2 = *it2;
                }

                // Control nodes must exist
                auto ctrlIt1 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == cn1; });
                if (ctrlIt1 == nodes.end()) {
                    throw Error("Error: Control node " + cn1 + " not found");
                }

                auto ctrlIt2 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == cn2; });
                if (ctrlIt2 == nodes.end()) {
                    throw Error("Error: Control node " + cn2 + " not found");
                }

                // Create and add VCCS
                double g = Coefficient_management(gain);
                VCCS* vccs = new VCCS(name, node2, node1, g, *ctrlIt1, *ctrlIt2);
                elements.push_back(vccs);
                success = true;
            } catch (const exception& e) {
                cout << e.what() << endl;
            }
        }

        void addCCVS(string name, string n1, string n2, string e, string gain) {
            try {
                for (auto el : elements) {
                    if (el->getName() == name) {
                        throw Error("Error: Element " + name + " already exists in the circuit");
                    }
                }

                // Find or create node1
                Node* node1 = new Node(n1, 0);
                auto it1 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n1; });
                if (it1 == nodes.end()) {
                    nodes.push_back(node1);
                } else {
                    delete node1;
                    node1 = *it1;
                }

                // Find or create node2
                Node* node2 = new Node(n2, 0);
                auto it2 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n2; });
                if (it2 == nodes.end()) {
                    nodes.push_back(node2);
                } else {
                    delete node2;
                    node2 = *it2;
                }

                // Find controlling element
                auto ctrlIt = std::find_if(elements.begin(), elements.end(), [&](Element* el) { return el->getName() == e; });
                if (ctrlIt == elements.end()) {
                    throw Error("Error: Controlling element " + e + " not found");
                }

                double g = Coefficient_management(gain);
                CCVS* ccvs = new CCVS(name, node2, node1, g, *ctrlIt);
                elements.push_back(ccvs);
                success = true;
            } catch (const exception& e) {
                cout << e.what() << endl;
            }
        }

        void addCCCS(string name, string n1, string n2, string e, string gain) {
            try {
                for (auto el : elements) {
                    if (el->getName() == name) {
                        throw Error("Error: Element " + name + " already exists in the circuit");
                    }
                }

                // Find or create node1
                Node* node1 = new Node(n1, 0);
                auto it1 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n1; });
                if (it1 == nodes.end()) {
                    nodes.push_back(node1);
                } else {
                    delete node1;
                    node1 = *it1;
                }

                // Find or create node2
                Node* node2 = new Node(n2, 0);
                auto it2 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n2; });
                if (it2 == nodes.end()) {
                    nodes.push_back(node2);
                } else {
                    delete node2;
                    node2 = *it2;
                }

                // Find controlling element
                auto ctrlIt = std::find_if(elements.begin(), elements.end(), [&](Element* el) { return el->getName() == e; });
                if (ctrlIt == elements.end()) {
                    throw Error("Error: Controlling element " + e + " not found");
                }

                double g = Coefficient_management(gain);
                CCCS* cccs = new CCCS(name, node2, node1, g, *ctrlIt);
                elements.push_back(cccs);
                success = true;
            } catch (const exception& e) {
                cout << e.what() << endl;
            }
        }

        void addDeltaV(string name, string n1, string n2, string tPulse, string area) {
            try {
                // Check for duplicate VoltageSource
                for (auto e : elements) {
                    if (e->getName() == name && dynamic_cast<VoltageSource*>(e)) {
                        throw Error("Error: VoltageSource " + name + " already exists in the circuit");
                    }
                }

                // Find or create node1
                Node* node1 = new Node(n1, 0);
                auto it1 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n1; });
                if (it1 == nodes.end()) {
                    nodes.push_back(node1);
                } else {
                    delete node1;
                    node1 = *it1;
                }

                // Find or create node2
                Node* node2 = new Node(n2, 0);
                auto it2 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n2; });
                if (it2 == nodes.end()) {
                    nodes.push_back(node2);
                } else {
                    delete node2;
                    node2 = *it2;
                }

                double value = unitHandler(area, "Voltage");
                double pulseTime = unitHandler2(tPulse);

                VoltageSource* Vs = new deltaVoltageSource(name, node2, node1, pulseTime, value);
                elements.push_back(Vs);
                success = true;
            } catch (const exception& e) {
                cout << e.what() << endl;
            }
        }

        void addDeltaI(string name, string n1, string n2, string tPulse, string area) {
            try {
                // Check for duplicate CurrentSource
                for (auto e : elements) {
                    if (e->getName() == name && dynamic_cast<CurrentSource*>(e)) {
                        throw Error("Error: CurrentSource " + name + " already exists in the circuit");
                    }
                }

                // Find or create node1
                Node* node1 = new Node(n1, 0);
                auto it1 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n1; });
                if (it1 == nodes.end()) {
                    nodes.push_back(node1);
                } else {
                    delete node1;
                    node1 = *it1;
                }

                // Find or create node2
                Node* node2 = new Node(n2, 0);
                auto it2 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n2; });
                if (it2 == nodes.end()) {
                    nodes.push_back(node2);
                } else {
                    delete node2;
                    node2 = *it2;
                }

                double value = unitHandler(area, "Current");
                double pulseTime = unitHandler2(tPulse);

                CurrentSource* Vs = new deltaCurrentSource(name, node2, node1, pulseTime, value);
                elements.push_back(Vs);
                success = true;
            } catch (const exception& e) {
                cout << e.what() << endl;
            }
        }

        void addSineVoltageSource(string name, string n1, string n2, string vOffset, string vAmplitude, string frequency) {
            try {
                for (auto e : elements) {
                    if (e->getName() == name && dynamic_cast<VoltageSource*>(e)) {
                        throw Error("Error: VoltageSource " + name + " already exists in the circuit");
                    }
                }

                Node* node1 = new Node(n1, 0);
                auto it1 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n1; });
                if (it1 == nodes.end()) {
                    nodes.push_back(node1);
                } else {
                    delete node1;
                    node1 = *it1;
                }

                Node* node2 = new Node(n2, 0);
                auto it2 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n2; });
                if (it2 == nodes.end()) {
                    nodes.push_back(node2);
                } else {
                    delete node2;
                    node2 = *it2;
                }

                double freq = unitHandler(frequency, "frequency");

                VoltageSource* vs = new SinVoltageSource(name, node2, node1, stod(vOffset), stod(vAmplitude), freq);
                elements.push_back(vs);
                success = true;
            } catch (const exception& e) {
                cout << e.what() << endl;
            }
        }

        void addSineCurrentSource(string name, string n1, string n2, string vOffset, string vAmplitude, string frequency) {
            try {
                for (auto e : elements) {
                    if (e->getName() == name && dynamic_cast<CurrentSource*>(e)) {
                        throw Error("Error: CurrentSource " + name + " already exists in the circuit");
                    }
                }

                Node* node1 = new Node(n1, 0);
                auto it1 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n1; });
                if (it1 == nodes.end()) {
                    nodes.push_back(node1);
                } else {
                    delete node1;
                    node1 = *it1;
                }

                Node* node2 = new Node(n2, 0);
                auto it2 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n2; });
                if (it2 == nodes.end()) {
                    nodes.push_back(node2);
                } else {
                    delete node2;
                    node2 = *it2;
                }

                double freq = unitHandler(frequency, "frequency");

                CurrentSource* cs = new SinCurrentSource(name, node2, node1, stod(vOffset), stod(vAmplitude), freq);
                elements.push_back(cs);
                success = true;
            } catch (const exception& e) {
                cout << e.what() << endl;
            }
        }

        void addPulseVoltageSource(string name, string n1, string n2, string Vinitial, string Von, string Tdelay, string Trise, string Tfall, string Ton, string Tperiod, string Ncycles = "0") {
            try {
                for (auto e : elements) {
                    if (e->getName() == name && dynamic_cast<VoltageSource*>(e)) {
                        throw Error("Error: VoltageSource " + name + " already exists in the circuit");
                    }
                }

                Node* node1 = new Node(n1, 0);
                auto it1 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n1; });
                if (it1 == nodes.end()) {
                    nodes.push_back(node1);
                } else {
                    delete node1;
                    node1 = *it1;
                }

                Node* node2 = new Node(n2, 0);
                auto it2 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n2; });
                if (it2 == nodes.end()) {
                    nodes.push_back(node2);
                } else {
                    delete node2;
                    node2 = *it2;
                }

                VoltageSource* vs = new PulseVoltageSource(name, node2, node1,
                                                           stod(Vinitial), stod(Von),
                                                           unitHandler2(Tdelay, "time"), unitHandler2(Trise, "time"),
                                                           unitHandler2(Tfall, "time"), unitHandler2(Ton, "time"),
                                                           unitHandler2(Tperiod, "time"), unitHandler2(Ncycles));

                elements.push_back(vs);
                success = true;
            } catch (const exception& e) {
                cout << e.what() << endl;
            }
        }

        void addPulseCurrentSource(string name, string n1, string n2, string Iinitial, string Ion, string Tdelay, string Trise, string Tfall, string Ton, string Tperiod, string Ncycles = "0") {
            try {
                for (auto e : elements) {
                    if (e->getName() == name && dynamic_cast<CurrentSource*>(e)) {
                        throw Error("Error: CurrentSource " + name + " already exists in the circuit");
                    }
                }

                Node* node1 = new Node(n1, 0);
                auto it1 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n1; });
                if (it1 == nodes.end()) {
                    nodes.push_back(node1);
                } else {
                    delete node1;
                    node1 = *it1;
                }

                Node* node2 = new Node(n2, 0);
                auto it2 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n2; });
                if (it2 == nodes.end()) {
                    nodes.push_back(node2);
                } else {
                    delete node2;
                    node2 = *it2;
                }

                CurrentSource* cs = new PulseCurrentSource(name, node2, node1,
                                                           stod(Iinitial), stod(Ion),
                                                           unitHandler2(Tdelay, "time"), unitHandler2(Trise, "time"),
                                                           unitHandler2(Tfall, "time"), unitHandler2(Ton, "time"),
                                                           unitHandler2(Tperiod, "time"), unitHandler2(Ncycles));

                elements.push_back(cs);
                success = true;
            } catch (const exception& e) {
                cout << e.what() << endl;
            }
        }

        void addDiode(string name, string n1, string n2, string model) {
            try {
                // Check for duplicate Diode
                for (auto e : elements) {
                    if (e->getName() == name && dynamic_cast<Diode*>(e)) {
                        throw Error("Error: Diode " + name + " already exists in the circuit");
                    }
                }

                // Find or create node1
                Node* node1 = new Node(n1, 0);
                auto it1 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n1; });
                if (it1 == nodes.end()) {
                    nodes.push_back(node1);
                } else {
                    delete node1;
                    node1 = *it1;
                }

                // Find or create node2
                Node* node2 = new Node(n2, 0);
                auto it2 = std::find_if(nodes.begin(), nodes.end(), [&](Node* n) { return n->getName() == n2; });
                if (it2 == nodes.end()) {
                    nodes.push_back(node2);
                } else {
                    delete node2;
                    node2 = *it2;
                }

                // Create diode based on model
                Diode* diode = nullptr;
                if (model == "Z") {
                    diode = new Diode(name, model, node2, node1, 0.7);
                } else if (model == "D") {
                    diode = new Diode(name, model, node2, node1, 0.0);
                } else {
                    throw Error("No diode with this model exists");
                }

                elements.push_back(diode);
                success = true;
            } catch (const exception& e) {
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

            // Initialize previous voltages and currents
            for (auto el : elements) {
                if (el->getType() == "Capacitor") {
                    capacitorVoltages[el->getName()] = 0.0;
                } else if (el->getType() == "Inductor") {
                    inductorCurrents[el->getName()] = 0.0;
                }
            }

            while (time <= Tstop + 1e-9) {
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

                        if (i != -1) G[i][i] += g;
                        if (j != -1) G[j][j] += g;
                        if (i != -1 && j != -1) {
                            G[i][j] -= g;
                            G[j][i] -= g;
                        }
                        if (i != -1) J[i] -= i_prev;
                        if (j != -1) J[j] += i_prev;
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

                vector<vector<double>> LU = LUDecompose(A);
                vector<double> X = LUSolve(LU, Z);

                TimePointResult result;
                result.time = time;
                for (auto it = nodeIndex.begin(); it != nodeIndex.end(); ++it)
                    result.nodeVoltages[it->first] = X[it->second];

                for (const auto& el : elements) {
                    if (el->getType() == "VoltageSource") {
                        int idx = voltageSourceIndex[el->getName()];
                        result.elementCurrents[el->getName()] = X[n + idx];
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
                        result.elementCurrents[name] = voltageAcross / el->getValue();
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

                for (auto el : elements) {
                    if (el->getType() == "Capacitor") {
                        int i = (el->getNode1()->getName() != "0" && el->getNode1()->getName() != "GND") ? nodeIndex[el->getNode1()->getName()] : -1;
                        int j = (el->getNode2()->getName() != "0" && el->getNode2()->getName() != "GND") ? nodeIndex[el->getNode2()->getName()] : -1;
                        double v1 = (i != -1) ? X[i] : 0;
                        double v2 = (j != -1) ? X[j] : 0;
                        capacitorVoltages[el->getName()] = v1 - v2;
                    } else if (el->getType() == "Inductor") {
                        int i = (el->getNode1()->getName() != "0" && el->getNode1()->getName() != "GND") ? nodeIndex[el->getNode1()->getName()] : -1;
                        int j = (el->getNode2()->getName() != "0" && el->getNode2()->getName() != "GND") ? nodeIndex[el->getNode2()->getName()] : -1;
                        double v1 = (i != -1) ? X[i] : 0;
                        double v2 = (j != -1) ? X[j] : 0;
                        double voltage = v1 - v2;
                        inductorCurrents[el->getName()] += (Tstep / el->getValue()) * voltage;
                    }
                }

                time += Tstep;
            }

            // Print
            vector<string> voltageNodes, currentElements;
            for (const auto& item : printItems) {
                if (item.size() > 2 && item.substr(0, 2) == "V(" && item.back() == ')') {
                    voltageNodes.push_back(item.substr(2, item.size() - 3));
                } else if (item.size() > 2 && item.substr(0, 2) == "I(" && item.back() == ')') {
                    currentElements.push_back(item.substr(2, item.size() - 3));
                }
            }

            for (const auto& n : voltageNodes) cout << "V(" << n << ")\t";
            for (const auto& e : currentElements) cout << "I(" << e << ")\t";
            cout << endl;

            for (const auto& res : simulationResults) {
                for (const auto& n : voltageNodes) cout << (res.nodeVoltages.count(n) ? to_string(res.nodeVoltages.at(n)) : "N/A") << "\t";
                for (const auto& e : currentElements) cout << (res.elementCurrents.count(e) ? to_string(res.elementCurrents.at(e)) : "N/A") << "\t";
                cout << endl;
            }
        }


        string analyzeTransientD(double tStart, double tEnd, double step, vector<Node*> wantedNodes, vector<Element*> wantedElements) {
            stringstream ss;
            map<string, int> nodeIndex;
            int idx = 0;

            for (Node* n : nodes) {
                if (!n->getIsGround()) {
                    nodeIndex[n->getName()] = idx++;
                }
            }
            int n = idx; // تعداد گره‌های غیر زمین
            int m = 0;   // تعداد منابع ولتاژ

            vector<VoltageSource*> voltageSources;
            for (Element* e : elements) {
                if (e->getType() == "VoltageSource") {
                    m++;
                    voltageSources.push_back(dynamic_cast<VoltageSource*>(e));
                }
            }

            int matrixSize = n + m;
            const int maxIter = 100;      // حداکثر تعداد تکرارها
            const double tolerance = 1e-6; // حد همگرایی

            // مقداردهی اولیه برای خازن و سلف
            for (Element* e : elements) {
                if (e->getType() == "Capacitor") {
                    Capacitor* c = dynamic_cast<Capacitor*>(e);
                    c->setPreviousVoltage(0.0);
                } else if (e->getType() == "Inductor") {
                    Inductor* i = dynamic_cast<Inductor*>(e);
                    i->setPreviousCurrent(0.0);
                }
            }

            // حل در هر گام زمانی
            for (double t = 0; t <= tEnd + 1e-12; t += step) {
                bool converged = false;
                vector<double> prevSolution;

                // حل تکرارشونده نیوتون-رافسون
                for (int iter = 0; iter < maxIter; ++iter) {
                    vector<vector<double>> G(n, vector<double>(n, 0.0));
                    vector<vector<double>> B(n, vector<double>(m, 0.0));
                    vector<vector<double>> D(m, vector<double>(n, 0.0));
                    vector<double> I(n, 0.0);
                    vector<double> F(m, 0.0);
                    int vsIndex = 0;

                    // بررسی هر المان
                    for (Element* e : elements) {
                        Node* pn = e->getNode1();
                        Node* nn = e->getNode2();
                        int pi = (!pn->getIsGround()) ? nodeIndex[pn->getName()] : -1;
                        int ni = (!nn->getIsGround()) ? nodeIndex[nn->getName()] : -1;

                        if (e->getType() == "Resistor") {
                            double r = e->getValue();
                            if (r <= 0) throw runtime_error("Resistor value must be positive");
                            double g = 1.0 / r;
                            if (pi != -1) G[pi][pi] += g;
                            if (ni != -1) G[ni][ni] += g;
                            if (pi != -1 && ni != -1) {
                                G[pi][ni] -= g;
                                G[ni][pi] -= g;
                            }
                        }
                        else if (e->getType() == "Capacitor") {
                            Capacitor* c = dynamic_cast<Capacitor*>(e);
                            double C = c->getValue();
                            if (C <= 0) throw runtime_error("Capacitor value must be positive");
                            double geq = C / step;
                            double prevV = c->getPreviousVoltage();
                            double ieq = geq * prevV;
                            if (pi != -1) G[pi][pi] += geq;
                            if (ni != -1) G[ni][ni] += geq;
                            if (pi != -1 && ni != -1) {
                                G[pi][ni] -= geq;
                                G[ni][pi] -= geq;
                            }
                            if (pi != -1) I[pi] += ieq;
                            if (ni != -1) I[ni] -= ieq;
                        }
                        else if (e->getType() == "Inductor") {
                            Inductor* i = dynamic_cast<Inductor*>(e);
                            double L = i->getValue();
                            if (L <= 0) throw runtime_error("Inductor value must be positive");
                            double geq = step / L;
                            double prevI = i->getPreviousCurrent();
                            if (pi != -1) G[pi][pi] += geq;
                            if (ni != -1) G[ni][ni] += geq;
                            if (pi != -1 && ni != -1) {
                                G[pi][ni] -= geq;
                                G[ni][pi] -= geq;
                            }
                            if (pi != -1) I[pi] -= prevI;
                            if (ni != -1) I[ni] += prevI;
                        }
                        else if (e->getType() == "VoltageSource") {
                            dynamic_cast<VoltageSource*>(e)->setValueAtTime(t);
                            if (pi != -1) B[pi][vsIndex] = 1;
                            if (ni != -1) B[ni][vsIndex] = -1;
                            if (pi != -1) D[vsIndex][pi] = 1;
                            if (ni != -1) D[vsIndex][ni] = -1;
                            F[vsIndex++] = e->getValue();
                        }
                        else if (e->getType() == "CurrentSource") {
                            dynamic_cast<CurrentSource*>(e)->setValueAtTime(t);
                            if (pi != -1) I[pi] += e->getValue();
                            if (ni != -1) I[ni] -= e->getValue();
                        }
                        else if (e->getType() == "DiodeD"||e->getType() == "DiodeZ") {
                            Diode* d = dynamic_cast<Diode*>(e);
                            d->updateState();

                            if (d->getIsOn()) {
                                // دیود روشن: معادل منبع ولتاژ با مقدار vOn
                                if (pi != -1) B[pi][vsIndex] = 1;
                                if (ni != -1) B[ni][vsIndex] = -1;
                                if (pi != -1) D[vsIndex][pi] = 1;
                                if (ni != -1) D[vsIndex][ni] = -1;
                                F[vsIndex++] = d->getValue();
                            }
                            // دیود خاموش: هیچ تاثیری ندارد (مدار باز)
                        }
                    }

                    // ساخت ماتریس کامل A و بردار Z
                    vector<vector<double>> A(matrixSize, vector<double>(matrixSize, 0.0));
                    vector<double> Z(matrixSize, 0.0);

                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < n; j++) A[i][j] = G[i][j];
                        for (int j = 0; j < m; j++) A[i][n + j] = B[i][j];
                        Z[i] = I[i];
                    }
                    for (int i = 0; i < m; i++) {
                        for (int j = 0; j < n; j++) A[n + i][j] = D[i][j];
                        Z[n + i] = F[i];
                    }

                    vector<vector<double>> LU = LUDecompose(A);
                    vector<double> solution = LUSolve(LU, Z);

                    // بررسی همگرایی
                    if (iter > 0) {
                        converged = true;
                        for (int i = 0; i < n; i++) {
                            if (abs(prevSolution[i] - solution[i]) > tolerance) {
                                converged = false;
                                break;
                            }
                        }
                        if (converged) break;
                    }

                    prevSolution = solution;

                    // به‌روزرسانی ولتاژ گره‌ها
                    for (Node* node : nodes) {
                        if (!node->getIsGround()) {
                            int nodeIdx = nodeIndex[node->getName()];
                            node->setVoltage(solution[nodeIdx]);
                        } else {
                            node->setVoltage(0.0);
                        }
                    }

                    // به‌روزرسانی جریان منابع ولتاژ
                    for (int i = 0; i < m; i++) {
                        voltageSources[i]->setCurrent(solution[n + i]);
                    }
                }

                if (!converged) {
                    throw runtime_error("Solution did not converge at time " + to_string(t));
                }

                // محاسبه جریان المان‌ها
                for (Element* e : elements) {
                    if (e->getType() == "Capacitor") {
                        Capacitor* c = dynamic_cast<Capacitor*>(e);
                        double v = c->getVoltage();
                        double prevV = c->getPreviousVoltage();
                        double current = (c->getValue() / step) * (v - prevV);
                        c->setCurrent(current);
                        c->setPreviousVoltage(v);
                    }
                    else if (e->getType() == "Inductor") {
                        Inductor* i = dynamic_cast<Inductor*>(e);
                        double v = i->getVoltage();
                        double current = i->getPreviousCurrent() + (step / i->getValue()) * v;
                        i->setCurrent(current);
                        i->setPreviousCurrent(current);
                    }
                    else if (e->getType() == "Resistor") {
                        e->setCurrent(e->getVoltage() / e->getValue());
                    }
                    else if (e->getType() == "Diode") {
                        Diode* d = dynamic_cast<Diode*>(e);
                        if (d->getIsOn()) {
                            // جریان دیود از حل سیستم به دست آمده است
                            // (به صورت جریان منبع ولتاژ متناظر)
                        } else {
                            d->setCurrent(0.0); // دیود خاموش: جریان صفر
                        }
                    }
                }

                // ذخیره نتایج در بازه مورد نظر
                if (t >= tStart) {
                    ss << fixed << setprecision(6);
                    ss << "Time = " << t << "s:\n";
                    for (auto n : wantedNodes) {
                        ss << "Node " << n->getName() << " voltage: " << n->getVoltage() << " V\n";
                    }
                    for (auto e : wantedElements) {
                        ss << e->getType() << " " << e->getName() << " current: " << e->getCurrent() << " A\n";
                    }
                }
            }
            return ss.str();
        }

        string tran(string tStart, string tEnd, string step,string wanted){
            double tStar =0;
            if (tStart != "0")
                tStar = unitHandler(tStart,"start time");
            double tStop = unitHandler(tEnd,"stop time");
            double ste = unitHandler(step,"time steps");

            vector<Node*>wantedNodes;
            vector<Element*>wantedElements;
            vector<string>words;

            stringstream ss(wanted);
            string word;
            while (ss>>word){
                if (word[0]=='V'){
                    Node* n=findNode(word.substr(2,word.length()-3));
                    if (!n){
                        throw Error("This node doesn't exist");
                    }
                    wantedNodes.push_back(findNode(word.substr(2,word.length()-3)));
                }
                else if (word[0]=='I'){
                    Element* e = findElement(word.substr(2,word.length()-3));
                    if (!e){
                        throw Error("This element doesn't exist");
                    }
                    wantedElements.push_back(findElement(word.substr(2,word.length()-3)));
                }
            }
            for (auto e:elements) {
                if (e->getType()=="DiodeD"||e->getType()=="DiodeZ"){
                    return analyzeTransientD(tStar,tStop,ste,wantedNodes,wantedElements);
                }
            }
            runTransient(ste, tStop, 0, -1, words);
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
                        cout << "2. NewFile <file_name>\n";
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
                            string fileName = command.substr(8);

                            ofstream newFile("schematics/" + fileName); // This will create the file if it doesn't exist
                            if (!newFile.is_open()) {
                                cout << "-Error: Cannot create file\n";
                                continue;
                            }


                            cout << "new file" << fileName << " created successfully" << endl;
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

//                        regex addResistorRegex(R"(\s*R (R\S+) (\S+) (\S+) (\S+)\s*)");
//                        regex addCapacitorRegex(R"(\s*C (C\S+) (\S+) (\S+) (\S+)\s*)");
//                        regex addInductorRegex(R"(\s*L (L\S+) (\S+) (\S+) (\S+)\s*)");
//                        regex addDiodeRegex(R"(\s*D (D\S+) (\S+) (\S+) (\S+)\s*)");
//                        regex addZDiodeRegex(R"(\s*Z (D\S+) (\S+) (\S+) (\S+)\s*)");
//                        regex addGroundRegex (R"(\s*GND (\S+)\s*)");
//                        regex addVCVSRegex (R"(E (E\S+) (\S+) (\S+) (\S+) (\S+) (\S+))");
//                        regex addVCCSRegex (R"(G (G\S+) (\S+) (\S+) (\S+) (\S+) (\S+))");
//                        regex addCCCSRegex (R"(F (F\S+) (\S+) (\S+) (\S+) (\S+))");
//                        regex addCCVSRegex (R"(H (H\S+) (\S+) (\S+) (\S+) (\S+))");
//                        regex addVoltageSourceRegex (R"(\s*VS (\S+) (\S+) (\S+) (\S+)\s*)");
//                        regex addSineVoltageSourceRegex (R"(\s*SINV (V\S+) (\S+) (\S+) SIN\((\S+) (\S+) (\S+)\)\s*)");
//                        regex addSineCurrentSourceRegex (R"(\s*SINI (I\S+) (\S+) (\S+) SIN\((\S+) (\S+) (\S+)\)\s*)");
//                        regex addPulseVoltageSourceRegex(R"(\s*PULSEV (V\S+) (\S+) (\S+) PULSE\s*\(\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)(?:\s+(\S+))?\s*\)\s*)");
//                        regex addPulseCurrentSourceRegex(R"(\s*PULSEI (I\S+) (\S+) (\S+) PULSE\s*\(\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)(?:\s+(\S+))?\s*\)\s*)");
//                        regex addCurrentSourceRegex (R"(\s*CS (\S+) (\S+) (\S+) (\S+)\s*)");
//                        regex printDCRegex(R"(\s*\.print DC (\S+) (\S+) (\S+) (\S+)(.*)\s*)");
//                        regex printOPRegex(R"(\s*\.print OP(.*)\s*)");
//                        regex printTranRegex(R"(\s*\.print TRAN (\S+) (\S+) (\S+)(.*)\s*)");
//                        smatch match;
//
//                        if (regex_match(input, match, addResistorRegex)) {
//                            circuit->addResistor(match[1], match[2], match[3], match[4]);
//                        }
//                        else if (regex_match(input,match,addVoltageSourceRegex)){
//                            circuit->addVoltageSource(match[1],match[2],match[3],match[4]);
//                        }
//                        else if (regex_match(input,match,addCurrentSourceRegex)){
//                            circuit->addCurrentSource(match[1],match[2],match[3],match[4]);
//                        }
//                        else if (regex_match(input, match, addCapacitorRegex)) {
//                            circuit->addCapacitor(match[1], match[2], match[3], match[4]);
//                        }
//                        else if (regex_match(input, match, addInductorRegex)) {
//                            circuit->addInductor(match[1], match[2], match[3], match[4]);
//                        }
//                        else if (regex_match(input, match, addDiodeRegex)) {
//                            circuit->addDiode(match[1], match[2], match[3], "D");
//                        }
//                        else if (regex_match(input, match, addZDiodeRegex)) {
//                            circuit->addDiode(match[1], match[2], match[3], "Z");
//                        }
//                        else if (regex_match(input, match, addVCVSRegex)){
//                            circuit->addVCVS(match[1],match[2],match[3],match[4],match[5],match[6]);
//                        }
//                        else if (regex_match(input, match, addCCVSRegex)){
//                            circuit->addCCVS(match[1],match[2],match[3],match[4],match[5]);
//                        }
//                        else if (regex_match(input, match, addVCCSRegex)){
//                            circuit->addVCCS(match[1],match[2],match[3],match[4],match[5],match[6]);
//                        }
//                        else if (regex_match(input, match, addCCCSRegex)){
//                            circuit->addCCCS(match[1],match[2],match[3],match[4],match[5]);
//                        }
//                        else if(regex_match(input, match, addPulseVoltageSourceRegex)) {
//                            if (match[11].matched) {
//                                circuit->addPulseVoltageSource(match[1], match[2], match[3],
//                                                      match[4], match[5], match[6],
//                                                      match[7], match[8], match[9],
//                                                      match[10], match[11]);
//                            } else {
//                                circuit->addPulseVoltageSource(match[1], match[2], match[3],
//                                                      match[4], match[5], match[6],
//                                                      match[7], match[8], match[9],
//                                                      match[10]);
//                            }
//                        }
//                        else if(regex_match(input, match, addPulseCurrentSourceRegex)) {
//                            if (match[11].matched) {
//                                circuit->addPulseCurrentSource(match[1], match[2], match[3],
//                                                      match[4], match[5], match[6],
//                                                      match[7], match[8], match[9],
//                                                      match[10], match[11]);
//                            } else {
//                                circuit->addPulseCurrentSource(match[1], match[2], match[3],
//                                                      match[4], match[5], match[6],
//                                                      match[7], match[8], match[9],
//                                                      match[10]);
//                            }
//                        }
//                        else if(regex_match(input,match,addSineVoltageSourceRegex)){
//                            circuit->addSineVoltageSource(match[1],match[2],match[3],match[4],match[5],match[6]);
//                        }
//                        else if(regex_match(input,match,addSineCurrentSourceRegex)){
//                            circuit->addSineCurrentSource(match[1],match[2],match[3],match[4],match[5],match[6]);
//                        }
//                        else if (regex_match(input, match, addPulseVoltageSourceRegex)) {
//                            circuit->addPulseVoltageSource(match[1], match[2], match[3], match[4], match[5],
//                                                  match[6], match[7], match[8], match[9], match[10]);
//                        }
//                        else if (regex_match(input, match, addPulseCurrentSourceRegex)) {
//                            circuit->addPulseCurrentSource(match[1], match[2], match[3], match[4], match[5],
//                                                  match[6], match[7], match[8], match[9], match[10]);
//                        }
//                        else if (regex_match(input,match,printTranRegex)){
//                            cout << circuit->tran(match[1],match[2],match[3],match[4]);
//                        } else {
//                            cout << "-Error: Syntax error\n";
//                        }
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

