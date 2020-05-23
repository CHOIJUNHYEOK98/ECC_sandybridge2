#pragma once
#include "ECC_CJH.h"

void Double_and_Add(IN word* k, IN point* P, IN word* p, OUT point* kP);
void generalMontgomery(IN word* k, IN point* P, IN word* p, OUT point* kP);
void windows(IN word* k, IN int w, IN word* p, OUT word* k2);
void precomputationWindow(IN point* P, IN int w, IN word* p, OUT point* pP);
void fixedWindow(IN word* k2, IN point* pP, IN int w, IN word* p, OUT point* kP);
