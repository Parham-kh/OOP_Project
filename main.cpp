#include <bits/stdc++.h>

using namespace std;
class ResistorInvalidAmount : public exception {
public:
    const char* what() const noexcept override {
        return "Error: Resistance cannot be zero or negative";
    }
};
class SyntaxError : public exception {
public:
    const char* what() const noexcept override {
        return "Error: Syntax error";
    }
};
namespace Model {
    class Node {
        private:
        string name="";
        double voltage=0;

        public:
        Node(string name, double voltage) :name(name),voltage(voltage){}
        string getName() {
            return name;
        }
    };
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
        virtual ~Element() = default;
    };
    class VoltageSource:public Element {};
    class Resistor:public Element {
        private:
        double resistance=0;

        public:
        Resistor(string name, Node* node1, Node* node2, double resistance)
        : Element(name, "Resistor", node1, node2, 0, 0), resistance(resistance) {}
        double getResistance() {return resistance;}
    };
    class Inductor:public Element{};
    class Capacitor:public Element{};
    class CurrentSource:public Element{};
};
using namespace Model;
namespace Controller {
    class Circuit {
        vector<Element*> elements;
        vector<Node*> nodes;

    public:
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
        void addResistor(string name,string n1,string n2,string resisrance) {
            try {
                Node* node1=new Node(n1,0);
                Node* node2=new Node(n2,0);
                nodes.push_back(node1);
                Resistor* resistor=new Resistor(name,node1,node2,Coefficient_management(resisrance));
                elements.push_back(resistor);
            }catch(const exception& e) {
                cout<<e.what()<<endl;
            }
        }
    };
}
using namespace Controller;
namespace View {
class view {
    private:
    Circuit* circuit;

    public:
    view(Circuit* c) : circuit(c) {}
    void run() {
        string input;
        getline(cin, input);
        regex r1("add R(\\S+) (\\S+) (\\S+) (\\S+)");
        smatch match;
        if (regex_search(input, match, r1)) {
            circuit->addResistor(match[1].str(), match[2].str(), match[3].str(), match[4].str());
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
