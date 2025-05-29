#include "View.h"

using namespace View;
using namespace std;

void view::run() {
    string input;
    getline(cin, input);
    regex r1("add R(\\S+) (\\S+) (\\S+) (\\S+)");
    regex r2("");
    smatch match;
    if (regex_search(input, match, r1)) {
        circuit->addResistor(match[1].str(), match[2].str(), match[3].str(), match[4].str());
    }
}
