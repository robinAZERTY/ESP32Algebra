/* this main file is used to test the compilation of the library */
#include <ekf.hpp>
#include <ukf.hpp>

void setup() {
}

void loop() {}


#ifdef NATIVE
int main()
{
    setup();
    while (true) loop();    
    return 0;
}
#endif