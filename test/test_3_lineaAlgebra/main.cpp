/*
to run all the test use the following command
pio test -e native
*/

#include <unity.h>
#include <linearAlgebra.hpp>

ldl_matrix<float> ldl_matrix_test;
triangMatrix<float> triang_matrix_test;

void setUp() {
    // Initialisation avant chaque test (laisser vide si inutile)
}

void tearDown() {
    // Nettoyage après chaque test (laisser vide si inutile)
}


void setup() {
    UNITY_BEGIN();
    UNITY_END();
}

void loop() {
    // Vide car Unity fonctionne avec setup() pour exécuter tous les tests une fois.
}


#ifdef NATIVE
int main(int argc, char **argv) {
    setup();
    return 0;
}
#endif