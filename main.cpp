#include <bits/stdc++.h>
using namespace std;
namespace Model {
    class Node{};
    class Element {
        string type;
    };
    class VoltageSource:public Element {
    };
    class Resistor:public Element {

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
    };
}
using namespace Controller;
namespace View {
    class view {
        public:
        void run() {
            string input;
            getline(cin, input);
            regex r1("");
        }
    };
}
using namespace View;


int main() {
    view Views;
    Views.run();
    return 0;
}