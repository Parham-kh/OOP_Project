#include <bits/stdc++.h>
#include <Exceptions/Exceptions.h>
#include <Model/Node.h>
#include <Model/Element.h>
#include <Model/VoltageSource.h>
#include <Model/Resistor.h>
#include <Model/Capacitor.h>
#include <Model/Inductor.h>
#include <Model/CurrentSource.h>
#include <ProgramController/Circuit.h>
#include <View/View.h>
using namespace std;
// class ResistorInvalidAmount : public exception {
// public:
//     const char* what() const noexcept override {
//         return "Error: Resistance cannot be zero or negative";
//     }
// };
// class SyntaxError : public exception {
// public:
//     const char* what() const noexcept override {
//         return "Error: Syntax error";
//     }
// };
// namespace Model {
//     // class Node {
//     //     private:
//     //     string name="";
//     //     double voltage=0;
//     //
//     //     public:
//     //     Node(string name, double voltage) :name(name),voltage(voltage){}
//     // };
//     // class Element {
//     //     private:
//     //     Node* node1;
//     //     Node* node2;
//     //
//     //     protected:
//     //     string name="";
//     //     string type="";
//     //     double voltage=0;
//     //     double current=0;
//     //
//     //     public:
//     //     Element(string name, string type, Node* node1, Node* node2, double voltage,double current)
//     //     : name(name), type(type), node1(node1), node2(node2), voltage(voltage), current(current) {}
//     //     virtual ~Element() = default;
//     // };
//     // class VoltageSource:public Element {};
//     // class Resistor:public Element {
//     //     private:
//     //     double resistance=0;
//     //
//     //     public:
//     //     Resistor(string name, Node* node1, Node* node2, double resistance)
//     //     : Element(name, "Resistor", node1, node2, 0, 0), resistance(resistance) {}
//     // };
//     // class Inductor:public Element{};
//     // class Capacitor:public Element{};
//     // class CurrentSource:public Element{};
// };
using namespace Model;
// namespace Controller {
//     // class Circuit {
//     //     vector<Element*> elements;
//     //     vector<Node*> nodes;
//     //
//     // public:
//     //     void addResistor(string name,string n1,string n2,double resisrance) {
//     //         try {
//     //             Node* node1=new Node(n1,0);
//     //             Node* node2=new Node(n2,0);
//     //             Resistor* resistor=new Resistor(name,node1,node2,resisrance);
//     //             elements.push_back(resistor);
//     //         }catch(const exception& e) {
//     //             cout<<e.what()<<endl;
//     //         }
//     //     }
//     // };
// }
using namespace Controller;
// namespace View {
// // class view {
// //     private:
// //     Circuit* circuit;
// //
// //     public:
// //     void run() {
// //         string input;
// //         getline(cin, input);
// //         regex r1("add R(\\S+) (\\S+) (\\S+) (\\S+)");
// //         smatch match;
// //         if (regex_search(input, match, r1)) {
// //             circuit->addResistor(match[1].str(), match[2].str(), match[3].str(), stod(match[4].str()));
// //         }
// //     }
// // };
// }
using namespace View;


int main() {
    Circuit circuit;
    view Views(&circuit);
    Views.run();
    return 0;
}