#include "View.h"

using namespace View;
using namespace std;

void view::run() {
    string input;
    getline(cin, input);
    regex r1("add R(\\S+) (\\S+) (\\S+) (\\S+)");
    smatch match;
    if (regex_search(input, match, r1)) {
        circuit->addResistor(match[1].str(), match[2].str(), match[3].str(), stod(match[4].str()));
    }
}
