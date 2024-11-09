/*
to run all the test use the following command
pio test -e native
*/

#include <unity.h>
#include <vector.hpp>
using namespace operators;
#ifdef NATIVE
template <typename T>
void print(const Vector<T> &v) {
    std::cout << to_string(v) << std::endl;
}
#else
template <typename T>
void print(const Vector<T> &v) {
    Serial.println(to_string(v).c_str());
}
#endif

void test_ctor_dtor(void) {
    Vector<double> v(3);
    v.fill(1.0);
    TEST_ASSERT_EQUAL(3, v.size());
    TEST_ASSERT_EQUAL(3, v.capacity());
    TEST_ASSERT_FALSE(v.shared());

    Vector<double> v2(v);
    TEST_ASSERT_EQUAL(3, v2.size());
    TEST_ASSERT_EQUAL(3, v2.capacity());
    TEST_ASSERT_FALSE(v2.shared());
    TEST_ASSERT_TRUE(v == v2);

    Vector<double> v3(v, true);
    TEST_ASSERT_EQUAL(3, v3.size());
    TEST_ASSERT_EQUAL(0, v3.capacity());
    TEST_ASSERT_TRUE(v3.shared());
    TEST_ASSERT_TRUE(v == v3);

    Vector<double> v4(&v3[0], 3);
    TEST_ASSERT_EQUAL(3, v4.size());
    TEST_ASSERT_EQUAL(3, v4.capacity());
    TEST_ASSERT_FALSE(v4.shared());
    TEST_ASSERT_TRUE(v == v4);

    TEST_ASSERT_EQUAL(0, internal::tmp<Vector<double>>::currentlyUsedCount());
}

void test_resize(void) {
    Vector<double> v(3);
    Vector<double> v2(3);
    v2 = v;

    v.resize(5);
    TEST_ASSERT_EQUAL(5, v.size());
    TEST_ASSERT_EQUAL(5, v.capacity());
    TEST_ASSERT_FALSE(v.shared());
    TEST_ASSERT_FALSE(v == v2);

    v.resize(3, false);
    TEST_ASSERT_EQUAL(3, v.size());
    TEST_ASSERT_EQUAL(5, v.capacity());
    TEST_ASSERT_FALSE(v.shared());
    TEST_ASSERT_FLOAT_WITHIN(1e-6, v[0], v2[0]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, v[1], v2[1]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, v[2], v2[2]);
}

void test_findFirst(void) {
    Vector<double> v(3);
    v[0] = 1;
    v[1] = 2;
    v[2] = 3;
    TEST_ASSERT_EQUAL(0, v.findFirst(1));
    TEST_ASSERT_EQUAL(1, v.findFirst(2));
    TEST_ASSERT_EQUAL(2, v.findFirst(3));
    TEST_ASSERT_EQUAL(-1, v.findFirst(4));
}

void test_findLast(void) {
    Vector<double> v(3);
    v[0] = 1;
    v[1] = 2;
    v[2] = 3;
    TEST_ASSERT_EQUAL(0, v.findLast(1));
    TEST_ASSERT_EQUAL(1, v.findLast(2));
    TEST_ASSERT_EQUAL(2, v.findLast(3));
    TEST_ASSERT_EQUAL(-1, v.findLast(4));
}

void test_findAll(void) {
    Vector<double> v(5);
    v[0] = 1;
    v[1] = 2;
    v[2] = 1;
    v[3] = 3;
    v[4] = 1;
    Vector<size_t> indices_holder;
    v.findAll(1, indices_holder);
    TEST_ASSERT_EQUAL(3, indices_holder.size());
    TEST_ASSERT_EQUAL(0, indices_holder[0]);
    TEST_ASSERT_EQUAL(2, indices_holder[1]);
    TEST_ASSERT_EQUAL(4, indices_holder[2]);
}

void test_removeFirst(void) {
    Vector<double> v(5);
    v[0] = 1;
    v[1] = 2;
    v[2] = 1;
    v[3] = 3;
    v[4] = 1;
    v.removeFirst(1);
    TEST_ASSERT_EQUAL(4, v.size());
    TEST_ASSERT_EQUAL(2, v[0]);
    TEST_ASSERT_EQUAL(1, v[1]);
    TEST_ASSERT_EQUAL(3, v[2]);
    TEST_ASSERT_EQUAL(1, v[3]);
}

void test_removeLast(void) {
    Vector<double> v(5);
    v[0] = 1;
    v[1] = 2;
    v[2] = 1;
    v[3] = 3;
    v[4] = 1;
    v.removeLast(1);
    TEST_ASSERT_EQUAL(4, v.size());
    TEST_ASSERT_EQUAL(1, v[0]);
    TEST_ASSERT_EQUAL(2, v[1]);
    TEST_ASSERT_EQUAL(1, v[2]);
    TEST_ASSERT_EQUAL(3, v[3]);
}

void test_removeAt(void) {
    Vector<double> v(5);
    v[0] = 1;
    v[1] = 2;
    v[2] = 1;
    v[3] = 3;
    v[4] = 1;
    v.removeAt(2);
    TEST_ASSERT_EQUAL(4, v.size());
    TEST_ASSERT_EQUAL(1, v[0]);
    TEST_ASSERT_EQUAL(2, v[1]);
    TEST_ASSERT_EQUAL(3, v[2]);
    TEST_ASSERT_EQUAL(1, v[3]);
}

void test_insert(void) {
    Vector<double> v(5);
    v[0] = 1;
    v[1] = 2;
    v[2] = 3;
    v[3] = 4;
    v[4] = 5;
    v.insert(6, 2);
    TEST_ASSERT_EQUAL(6, v.size());
    TEST_ASSERT_EQUAL(1, v[0]);
    TEST_ASSERT_EQUAL(2, v[1]);
    TEST_ASSERT_EQUAL(6, v[2]);
    TEST_ASSERT_EQUAL(3, v[3]);
    TEST_ASSERT_EQUAL(4, v[4]);
    TEST_ASSERT_EQUAL(5, v[5]);
}

void test_push_back(void) {
    Vector<double> v(5);
    v[0] = 1;
    v[1] = 2;
    v[2] = 3;
    v[3] = 4;
    v[4] = 5;
    v.push_back(6);
    TEST_ASSERT_EQUAL(6, v.size());
    TEST_ASSERT_EQUAL(1, v[0]);
    TEST_ASSERT_EQUAL(2, v[1]);
    TEST_ASSERT_EQUAL(3, v[2]);
    TEST_ASSERT_EQUAL(4, v[3]);
    TEST_ASSERT_EQUAL(5, v[4]);
    TEST_ASSERT_EQUAL(6, v[5]);
}

void test_pop_back(void) {
    Vector<double> v(5);
    v[0] = 1;
    v[1] = 2;
    v[2] = 3;
    v[3] = 4;
    v[4] = 5;
    v.pop_back();
    TEST_ASSERT_EQUAL(4, v.size());
    TEST_ASSERT_EQUAL(1, v[0]);
    TEST_ASSERT_EQUAL(2, v[1]);
    TEST_ASSERT_EQUAL(3, v[2]);
    TEST_ASSERT_EQUAL(4, v[3]);
}

void test_sort(void) {
    Vector<double> v(5);
    v[0] = 5;
    v[1] = 4;
    v[2] = 3;
    v[3] = 2;
    v[4] = 1;
    v.sort();
    TEST_ASSERT_EQUAL(5, v.size());
    TEST_ASSERT_EQUAL(1, v[0]);
    TEST_ASSERT_EQUAL(2, v[1]);
    TEST_ASSERT_EQUAL(3, v[2]);
    TEST_ASSERT_EQUAL(4, v[3]);
    TEST_ASSERT_EQUAL(5, v[4]);
}

void test_fill(void) {
    Vector<double> v(5);
    v.fill(1);
    TEST_ASSERT_EQUAL(5, v.size());
    TEST_ASSERT_EQUAL(1, v[0]);
    TEST_ASSERT_EQUAL(1, v[1]);
    TEST_ASSERT_EQUAL(1, v[2]);
    TEST_ASSERT_EQUAL(1, v[3]);
    TEST_ASSERT_EQUAL(1, v[4]);
}

void test_hold(void) {
    Vector<double> v(5);
    v.fill(1);
    Vector<double> v2(5);
    v2.hold(v);
    TEST_ASSERT_EQUAL(5, v2.size());
    TEST_ASSERT_EQUAL(1, v2[0]);
    TEST_ASSERT_EQUAL(1, v2[1]);
    TEST_ASSERT_EQUAL(1, v2[2]);
    TEST_ASSERT_EQUAL(1, v2[3]);
    TEST_ASSERT_EQUAL(1, v2[4]);
}

void test_holdAdd(void) {
    Vector<double> v(5);
    v.fill(1);
    Vector<double> v2(5);
    v2.fill(2);
    Vector<double> v3(5);
    v3.holdAdd(v, v2);
    TEST_ASSERT_EQUAL(5, v3.size());
    TEST_ASSERT_EQUAL(3, v3[0]);
    TEST_ASSERT_EQUAL(3, v3[1]);
    TEST_ASSERT_EQUAL(3, v3[2]);
    TEST_ASSERT_EQUAL(3, v3[3]);
    TEST_ASSERT_EQUAL(3, v3[4]);
}

void test_holdSub(void) {
    Vector<double> v(5);
    v.fill(1);
    Vector<double> v2(5);
    v2.fill(2);
    Vector<double> v3(5);
    v3.holdSub(v, v2);
    TEST_ASSERT_EQUAL(5, v3.size());
    TEST_ASSERT_EQUAL(-1, v3[0]);
    TEST_ASSERT_EQUAL(-1, v3[1]);
    TEST_ASSERT_EQUAL(-1, v3[2]);
    TEST_ASSERT_EQUAL(-1, v3[3]);
    TEST_ASSERT_EQUAL(-1, v3[4]);
}

void test_holdMul(void) {
    Vector<double> v(5);
    v.fill(1);
    Vector<double> v2(5);
    v2.fill(2);
    Vector<double> v3(5);
    v3.holdMul(v, v2);
    TEST_ASSERT_EQUAL(5, v3.size());
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 2, v3[0]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 2, v3[1]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 2, v3[2]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 2, v3[3]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 2, v3[4]);
}

void test_holdDiv(void) {
    Vector<double> v(5);
    v.fill(1);
    Vector<double> v2(5);
    v2.fill(2);
    Vector<double> v3(5);
    v3.holdDiv(v, v2);
    TEST_ASSERT_EQUAL(5, v3.size());
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 0.5, v3[0]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 0.5, v3[1]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 0.5, v3[2]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 0.5, v3[3]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 0.5, v3[4]);
}

void test_holdAdd2(void) {
    Vector<double> v(5);
    v.fill(1);
    Vector<double> v2(5);
    v2.fill(2);
    Vector<double> v3(5);
    v3.holdAdd(v, 2);
    TEST_ASSERT_EQUAL(5, v3.size());
    TEST_ASSERT_EQUAL(3, v3[0]);
    TEST_ASSERT_EQUAL(3, v3[1]);
    TEST_ASSERT_EQUAL(3, v3[2]);
    TEST_ASSERT_EQUAL(3, v3[3]);
    TEST_ASSERT_EQUAL(3, v3[4]);
}

void test_holdSub2(void) {
    Vector<double> v(5);
    v.fill(1);
    Vector<double> v2(5);
    v2.fill(2);
    Vector<double> v3(5);
    v3.holdSub(v, 2);
    TEST_ASSERT_EQUAL(5, v3.size());
    TEST_ASSERT_EQUAL(-1, v3[0]);
    TEST_ASSERT_EQUAL(-1, v3[1]);
    TEST_ASSERT_EQUAL(-1, v3[2]);
    TEST_ASSERT_EQUAL(-1, v3[3]);
    TEST_ASSERT_EQUAL(-1, v3[4]);
}

void test_holdMul2(void) {
    Vector<double> v(5);
    v.fill(1);
    Vector<double> v3(5);
    v3.holdMul(v, 2);
    TEST_ASSERT_EQUAL(5, v3.size());
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 2, v3[0]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 2, v3[1]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 2, v3[2]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 2, v3[3]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 2, v3[4]);
}

void test_holdDiv2(void) {
    Vector<double> v(5);
    v.fill(1);
    Vector<double> v3(5);
    // v3.holdDiv(v, 2.0);
    // TEST_ASSERT_EQUAL(5, v3.size());
    // TEST_ASSERT_FLOAT_WITHIN(1e-6, 0.5, v3[0]);
    // TEST_ASSERT_FLOAT_WITHIN(1e-6, 0.5, v3[1]);
    // TEST_ASSERT_FLOAT_WITHIN(1e-6, 0.5, v3[2]);
    // TEST_ASSERT_FLOAT_WITHIN(1e-6, 0.5, v3[3]);
    // TEST_ASSERT_FLOAT_WITHIN(1e-6, 0.5, v3[4]);
}

void test_holdSub3(void) {
    Vector<double> v(5);
    v.fill(1);
    Vector<double> v3(5);
    v3.holdSub(2, v);
    TEST_ASSERT_EQUAL(5, v3.size());
    TEST_ASSERT_EQUAL(1, v3[0]);
    TEST_ASSERT_EQUAL(1, v3[1]);
    TEST_ASSERT_EQUAL(1, v3[2]);
    TEST_ASSERT_EQUAL(1, v3[3]);
    TEST_ASSERT_EQUAL(1, v3[4]);
}

void test_holdDiv3(void) {
    Vector<double> v(5);
    v.fill(1);
    Vector<double> v3(5);
    v3.holdDiv(2, v);
    TEST_ASSERT_EQUAL(5, v3.size());
    TEST_ASSERT_EQUAL(2, v3[0]);
    TEST_ASSERT_EQUAL(2, v3[1]);
    TEST_ASSERT_EQUAL(2, v3[2]);
    TEST_ASSERT_EQUAL(2, v3[3]);
    TEST_ASSERT_EQUAL(2, v3[4]);
}

void test_equal(void) {
    Vector<double> v(5);
    v.fill(1);
    Vector<double> v2(5);
    v2.fill(1);
    TEST_ASSERT_TRUE(v == v2);
}

void test_notEqual(void) {
    Vector<double> v(5);
    v.fill(1);
    Vector<double> v2(5);
    v2.fill(2);
    TEST_ASSERT_TRUE(v != v2);
}

void test_equalAffect(void) {
    Vector<double> v(5);
    v.fill(1);
    Vector<double> v2;
    v2 = v;
    TEST_ASSERT_TRUE(v == v2);
}

void test_addEqual(void) {
    Vector<double> v(5);
    v.fill(1);
    Vector<double> v2(5);
    v2.fill(2);
    v += v2;
    TEST_ASSERT_EQUAL(5, v.size());
    TEST_ASSERT_EQUAL(3, v[0]);
    TEST_ASSERT_EQUAL(3, v[1]);
    TEST_ASSERT_EQUAL(3, v[2]);
    TEST_ASSERT_EQUAL(3, v[3]);
    TEST_ASSERT_EQUAL(3, v[4]);
}

void test_subEqual(void) {
    Vector<double> v(5);
    v.fill(1);
    Vector<double> v2(5);
    v2.fill(2);
    v -= v2;
    TEST_ASSERT_EQUAL(5, v.size());
    TEST_ASSERT_EQUAL(-1, v[0]);
    TEST_ASSERT_EQUAL(-1, v[1]);
    TEST_ASSERT_EQUAL(-1, v[2]);
    TEST_ASSERT_EQUAL(-1, v[3]);
    TEST_ASSERT_EQUAL(-1, v[4]);
}

void test_mulEqual(void) {
    Vector<double> v(5);
    v.fill(1);
    Vector<double> v2(5);
    v2.fill(2);
    v *= v2;
    TEST_ASSERT_EQUAL(5, v.size());
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 2, v[0]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 2, v[1]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 2, v[2]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 2, v[3]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 2, v[4]);
}

void test_divEqual(void) {
    Vector<double> v(5);
    v.fill(1);
    Vector<double> v2(5);
    v2.fill(2);
    v /= v2;
    TEST_ASSERT_EQUAL(5, v.size());
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 0.5, v[0]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 0.5, v[1]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 0.5, v[2]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 0.5, v[3]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 0.5, v[4]);
}

void test_equal2(void) {
    Vector<double> v(5);
    v.fill(1);
    TEST_ASSERT_TRUE(v == 1);
}

void test_notEqual2(void) {
    Vector<double> v(5);
    v.fill(1);
    TEST_ASSERT_TRUE(v != 2);
}

void test_equalAffect2(void) {
    Vector<double> v(5);
    v.fill(1);
    Vector<double> v2(5);
    v2.fill(1);
    TEST_ASSERT_TRUE(v == v2);
}

void test_addEqual2(void) {
    Vector<double> v(5);
    v.fill(1);
    v += 2;
    TEST_ASSERT_EQUAL(5, v.size());
    TEST_ASSERT_EQUAL(3, v[0]);
    TEST_ASSERT_EQUAL(3, v[1]);
    TEST_ASSERT_EQUAL(3, v[2]);
    TEST_ASSERT_EQUAL(3, v[3]);
    TEST_ASSERT_EQUAL(3, v[4]);
}

void test_subEqual2(void) {
    Vector<double> v(5);
    v.fill(1);
    v -= 2;
    TEST_ASSERT_EQUAL(5, v.size());
    TEST_ASSERT_EQUAL(-1, v[0]);
    TEST_ASSERT_EQUAL(-1, v[1]);
    TEST_ASSERT_EQUAL(-1, v[2]);
    TEST_ASSERT_EQUAL(-1, v[3]);
    TEST_ASSERT_EQUAL(-1, v[4]);
}

void test_mulEqual2(void) {
    Vector<double> v(5);
    v.fill(1);
    v *= 2;
    TEST_ASSERT_EQUAL(5, v.size());
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 2, v[0]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 2, v[1]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 2, v[2]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 2, v[3]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 2, v[4]);
}

void test_divEqual2(void) {
    Vector<double> v(5);
    v.fill(1);
    v /= 2;
    TEST_ASSERT_EQUAL(5, v.size());
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 0.5, v[0]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 0.5, v[1]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 0.5, v[2]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 0.5, v[3]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 0.5, v[4]);
}

void test_addoperator(void) {
    Vector<float> v(5);
    v.fill(1);
    v[0] += 2;
    Vector<float> v2(5);
    v2.fill(2);
    Vector<float> v3 = v + v2;
    TEST_ASSERT_EQUAL(5, v3.size());
    TEST_ASSERT_EQUAL(5, v3[0]);
    TEST_ASSERT_EQUAL(3, v3[1]);
    TEST_ASSERT_EQUAL(3, v3[2]);
    TEST_ASSERT_EQUAL(3, v3[3]);
    TEST_ASSERT_EQUAL(3, v3[4]);

    v3 = v + 2.0f;
    TEST_ASSERT_EQUAL(5, v3.size());
    TEST_ASSERT_EQUAL(5, v3[0]);
    TEST_ASSERT_EQUAL(3, v3[1]);
    TEST_ASSERT_EQUAL(3, v3[2]);
    TEST_ASSERT_EQUAL(3, v3[3]);
    TEST_ASSERT_EQUAL(3, v3[4]);
    
    Vector<float> v4 = 2.0f + v;
    v4 += 0.5f + v2 + 0.5f;
    TEST_ASSERT_EQUAL(5, v4.size());
    TEST_ASSERT_EQUAL(8, v4[0]);
    TEST_ASSERT_EQUAL(6, v4[1]);
    TEST_ASSERT_EQUAL(6, v4[2]);
    TEST_ASSERT_EQUAL(6, v4[3]);
    TEST_ASSERT_EQUAL(6, v4[4]);
}

void test_suboperator(void) {
    Vector<double> v(5);
    v.fill(1);
    v[0] -= 2;
    Vector<double> v2(5);
    v2.fill(2);
    Vector<double> v3 = v - v2;
    TEST_ASSERT_EQUAL(5, v3.size());
    TEST_ASSERT_EQUAL(-3, v3[0]);
    TEST_ASSERT_EQUAL(-1, v3[1]);
    TEST_ASSERT_EQUAL(-1, v3[2]);
    TEST_ASSERT_EQUAL(-1, v3[3]);
    TEST_ASSERT_EQUAL(-1, v3[4]);
}

void test_muloperator(void) {
    Vector<float> v(5);
    v.fill(1);
    v[0] *= 2;
    Vector<float> v2(5);
    v2.fill(2);
    Vector<float> v3(5);
    v3.holdMul(v, v2);
    TEST_ASSERT_EQUAL(5, v3.size());
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 4, v3[0]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 2, v3[1]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 2, v3[2]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 2, v3[3]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 2, v3[4]);
    
    v3 = v * 2.0f;
    TEST_ASSERT_EQUAL(5, v3.size());
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 4, v3[0]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 2, v3[1]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 2, v3[2]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 2, v3[3]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 2, v3[4]);
    
    Vector<float> v4 = 2.0f * v;
    v4 *= 0.5f * v2 * 0.5f;
    TEST_ASSERT_EQUAL(5, v4.size());
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 2, v4[0]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 1, v4[1]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 1, v4[2]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 1, v4[3]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 1, v4[4]);
}

void test_divoperator(void) {
    Vector<double> v(5);
    v.fill(1);
    v[0] /= 2;
    Vector<double> v2(5);
    v2.fill(2);
    Vector<double> v3 = v / v2;
    TEST_ASSERT_EQUAL(5, v3.size());
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 0.25, v3[0]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 0.5, v3[1]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 0.5, v3[2]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 0.5, v3[3]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 0.5, v3[4]);

    v3 = v / 2.0;
    TEST_ASSERT_EQUAL(5, v3.size());
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 0.25, v3[0]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 0.5, v3[1]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 0.5, v3[2]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 0.5, v3[3]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 0.5, v3[4]);
    
    Vector<double> v4 = 2.0 / v;
    v4 /= 0.5 / v2 / 0.5;
    TEST_ASSERT_EQUAL(5, v4.size());
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 8, v4[0]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 4, v4[1]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 4, v4[2]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 4, v4[3]);
    TEST_ASSERT_FLOAT_WITHIN(1e-6, 4, v4[4]);
}



// Ajoute des fonctions de test supplémentaires ici, en utilisant le même format.

void setUp() {
    // Initialisation avant chaque test (laisser vide si inutile)
}

void tearDown() {
    // Nettoyage après chaque test (laisser vide si inutile)
}


void setup() {
    UNITY_BEGIN();
    RUN_TEST(test_ctor_dtor);
    RUN_TEST(test_resize);
    RUN_TEST(test_findFirst);
    RUN_TEST(test_findLast);
    RUN_TEST(test_findAll);
    RUN_TEST(test_removeFirst);
    RUN_TEST(test_removeLast);
    RUN_TEST(test_removeAt);
    RUN_TEST(test_insert);
    RUN_TEST(test_push_back);
    RUN_TEST(test_pop_back);
    RUN_TEST(test_sort);
    RUN_TEST(test_fill);
    RUN_TEST(test_hold);
    RUN_TEST(test_holdAdd);
    RUN_TEST(test_holdSub);
    RUN_TEST(test_holdMul);
    RUN_TEST(test_holdDiv);
    RUN_TEST(test_holdAdd2);
    RUN_TEST(test_holdSub2);
    RUN_TEST(test_holdMul2);
    RUN_TEST(test_holdDiv2);
    RUN_TEST(test_holdSub3);
    RUN_TEST(test_holdDiv3);
    RUN_TEST(test_equal);
    RUN_TEST(test_notEqual);
    RUN_TEST(test_equalAffect);
    RUN_TEST(test_addEqual);
    RUN_TEST(test_subEqual);
    RUN_TEST(test_mulEqual);
    RUN_TEST(test_divEqual);
    RUN_TEST(test_equal2);
    RUN_TEST(test_notEqual2);
    RUN_TEST(test_equalAffect2);
    RUN_TEST(test_addEqual2);
    RUN_TEST(test_subEqual2);
    RUN_TEST(test_mulEqual2);
    RUN_TEST(test_divEqual2);
    RUN_TEST(test_addoperator);
    RUN_TEST(test_suboperator);
    RUN_TEST(test_muloperator);
    RUN_TEST(test_divoperator);

    // Ajoute d'autres tests ici.
    #ifdef NATIVE
    UNITY_END();
    #endif
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