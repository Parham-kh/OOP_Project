#include <bits/stdc++.h>
using namespace std;

class Error : public exception {
    string message;

public:
    Error(const string& msg) : message(msg) {}
    const char* what() const noexcept override {
        return message.c_str();
    }
};

namespace Model {
    class Node {
        private:
        string name="";
        double voltage=0;
        bool IS_GND=false;

        public:
        Node(string name, double voltage) :name(name),voltage(voltage){}
        string getName() {
            return name;
        }
        void setGND(bool IS_GND) {
            this->IS_GND = IS_GND;
        }
        bool isGND() {
            return this->IS_GND;
        }
    };

    class Element {
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
        virtual ~Element() = default;
    };

    class VoltageSource:public Element {};

    class Resistor:public Element {
        double resistance = 0;

    public:
        Resistor(string name, Node* node1, Node* node2, double resistance)
        : Element(name, "Resistor", node1, node2, 0, 0), resistance(resistance) {}
    };

    class Inductor:public Element{
        double inductance = 0;

    public:
        Inductor(string name, Node* node1, Node* node2, double inductance)
        : Element(name, "Resistor", node1, node2, 0, 0), inductance(inductance) {}
    };

    class Capacitor : public Element {
        double capacitance = 0;

    public:
        Capacitor(string name, Node* node1, Node* node2, double capacitance)
                : Element(name, "Resistor", node1, node2, 0, 0), capacitance(capacitance) {}
        double getCapacitance() {
            return capacitance;
        }
    };

    class CurrentSource:public Element{};
    class Diode:public Element {
        string model="";

    public:
        Diode(string name, Node* node1, Node* node2, string model)
        : Element(name, "Diode", node1, node2, 0, 0), model(model) {}
    };
};
using namespace Model;

namespace Controller {
    class Circuit {
        vector<Element*> elements;
        vector<Node*> nodes;

    public:
        double Coefficient_management(string a) {
            try {
                regex r1("(\\S*)([k|m|u|n|Meg])");
                smatch match;
                if(regex_search(a, match, r1)) {
                    double output = stod(match[1].str());
                    if(output <= 0)
                        throw Error("Error: Resistance cannot be zero or negative");
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
                    if(stod(a) <= 0)
                        throw Error("Error: Resistance cannot be zero or negative");
                    return stod(a);
                }
            }
            catch(const exception& e) {
                cout<<e.what()<<endl;
            }
        }

        void addResistor(string name,string n1,string n2,string resisrance) {
            try {
                for(auto e:elements) {
                    if(e->getName()==name) {
                        throw Error("Error: Resistor "+name+" already exists in the circuit");
                    }
                }
                Node* node1;
                Node* node2;
                bool nod1=false,nod2=false;
                for(auto n:nodes) {
                    if(n->getName()==n1) {
                        nod1=true;
                        node1=n;
                    }
                    if(n->getName()==n2) {
                        nod2=true;
                        node2=n;
                    }
                }
                if(!nod1) {
                    Node* node1=new Node(n1,0);
                    nodes.push_back(node1);
                }
                if(!nod2) {
                    Node* node2=new Node(n2,0);
                    nodes.push_back(node2);
                }
                Resistor* resistor=new Resistor(name,node1,node2,Coefficient_management(resisrance));
                elements.push_back(resistor);
            }catch(const exception& e) {
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
                Node* node1;
                Node* node2;
                bool nod1=false,nod2=false;
                for(auto n:nodes) {
                    if(n->getName()==n1) {
                        nod1=true;
                        node1=n;
                    }
                    if(n->getName()==n2) {
                        nod2=true;
                        node2=n;
                    }
                }
                if(!nod1) {
                    Node* node1=new Node(n1,0);
                    nodes.push_back(node1);
                }
                if(!nod2) {
                    Node* node2=new Node(n2,0);
                    nodes.push_back(node2);
                }
                Capacitor* capacitor = new Capacitor(name, node1, node2, Coefficient_management(capacitance));
                elements.push_back(capacitor);
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
                Node* node1;
                Node* node2;
                bool nod1=false,nod2=false;
                for(auto n:nodes) {
                    if(n->getName()==n1) {
                        nod1=true;
                        node1=n;
                    }
                    if(n->getName()==n2) {
                        nod2=true;
                        node2=n;
                    }
                }
                if(!nod1) {
                    Node* node1=new Node(n1,0);
                    nodes.push_back(node1);
                }
                if(!nod2) {
                    Node* node2=new Node(n2,0);
                    nodes.push_back(node2);
                }
                Inductor* inductor = new Inductor(name, node1, node2, Coefficient_management(inductance));
                elements.push_back(inductor);
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

        void addDiode(string name,string n1,string n2,string model) {
            try {
                for(auto e:elements) {
                    if(e->getName()==name) {
                        throw Error("Error: Diode "+name+" already exists in the circuit");
                    }
                }
                Node* node1;
                Node* node2;
                bool nod1=false,nod2=false;
                for(auto n:nodes) {
                    if(n->getName()==n1) {
                        nod1=true;
                        node1=n;
                    }
                    if(n->getName()==n2) {
                        nod2=true;
                        node2=n;
                    }
                }
                if(!nod1) {
                    Node* node1=new Node(n1,0);
                    nodes.push_back(node1);
                }
                if(!nod2) {
                    Node* node2=new Node(n2,0);
                    nodes.push_back(node2);
                }
                Diode* diode=new Diode(name,node1,node2,model);
                elements.push_back(diode);
            }catch(const exception& e) {
                cout<<e.what()<<endl;
            }
        }

        void deleteDiode(string name) {
            bool found=false;
            for(int i=0;i<elements.size();i++) {
                if(elements[i]->getName()==name) {
                    elements.erase(elements.begin()+i);
                    found=true;
                }
            }
            if(!found) {
                throw Error("Error: Cannot delete diode; component not found");
            }
        }

        void add_GND(string node) {
            bool found=false;
            for(auto n:nodes) {
                if(n->getName()==node) {
                    n->setGND(true);
                    found=true;
                }
            }
            if(!found) {
                cout<<"Node does not exist"<<endl;
            }
        }

        void delete_GND(string node) {
            bool found=false;
            for(auto n:nodes) {
                if(n->getName()==node) {
                    n->setGND(false);
                    found=true;
                }
            }
            if(!found) {
                cout<<"Node does not exist"<<endl;
            }
        }
    };
}

using namespace Controller;

namespace View {
    class view {
    private:
        Circuit *circuit;

    public:
        view(Circuit *c) : circuit(c) {}

        void run() {
            while(true){
                try {
                    string input;
                    getline(cin, input);
                    regex add("add (\\S+) (\\S+) (\\S+) (\\S+)");
                    regex del("delete (\\S+)");
                    regex add_GND("add (\\S+) (\\S+)");
                    regex Delete_GND("delete (\\S+) (\\S+)");

                    smatch match;
                    if (regex_search(input, match, add)) {
                        if (match[1].str()[0] == 'R')
                            circuit->addResistor(match[1].str(), match[2].str(), match[3].str(), match[4].str());
                        else if (match[1].str()[0] == 'C')
                            circuit->addCapacitor(match[1].str(), match[2].str(), match[3].str(), match[4].str());
                        else if (match[1].str()[0] == 'L')
                            circuit->addInductor(match[1].str(), match[2].str(), match[3].str(), match[4].str());
                        else if (match[1].str()[0] == 'D')
                            circuit->addDiode(match[1].str(), match[2].str(), match[3].str(), match[4].str());
                        else
                            throw Error("Error: Element "+ match[1].str() + " not found in library");
                    }
                    else if (regex_search(input, match, del)) {
                        if (match[1].str()[0] == 'R')
                            circuit->deleteResistor(match[1].str());
                        else if (match[1].str()[0] == 'C')
                            circuit->deleteCapacitor(match[1].str());
                        else if (match[1].str()[0] == 'L')
                            circuit->deleteInductor(match[1].str());
                        else if (match[1].str()[0] == 'D')
                            circuit->deleteDiode(match[1].str());
                        else
                            throw Error("Error: Element "+ match[1].str() + " not found in library");
                    }
                    else if(regex_search(input, match, add_GND)) {
                        if(match[1].str() == "GND")
                            circuit->add_GND(match[2].str());
                        else
                            throw Error("Error: Element "+match[1].str() + " not found in library");
                    }
                    else if(regex_search(input, match, Delete_GND)) {
                        if(match[1].str() == "GND")
                            circuit->delete_GND(match[2].str());
                        else
                            throw Error("Error: Element "+match[1].str() + " not found in library");
                    }
                    else {
                        throw Error("Error: Syntax error");
                    }
                }
                catch(const exception& e) {
                    cout<<e.what()<<endl;
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