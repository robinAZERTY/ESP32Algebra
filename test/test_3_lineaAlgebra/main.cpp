/*
to run all the test use the following command
pio test -e native
*/

#include <unity.h>
#include <linearAlgebra.hpp>

void setUp() {
    // Initialisation avant chaque test (laisser vide si inutile)
}

void tearDown() {
    // Nettoyage après chaque test (laisser vide si inutile)
}


void setup() {
    UNITY_BEGIN();
    // Ajoute d'autres tests ici.
    UNITY_END();
}

void loop() {
    // Vide car Unity fonctionne avec setup() pour exécuter tous les tests une fois.
}

#ifdef NATIVE
int main(int argc, char **argv) {
    setup();
    return UNITY_END();
}
#endif