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
#include <vector>

using namespace std;
using namespace Model;
using namespace Controller;
using namespace View;


int main() {
    view Views;
    Views.run();
    return 0;
}