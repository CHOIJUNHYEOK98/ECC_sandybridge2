#include "ECC_CJH.h"

void Double_and_Add(IN word* k, IN point* P, IN word* p, OUT point* kP)
{
	int i, j;
	point R;
	R.flag = 1;

	for (i = SIZE-1; i >= 0; i--)
	{
		for (j = 31; j >= 0; j--)
		{
			R.flag = jacobianECDBL(&R, p, &R);
			if ((k[i] >> j) & 0x1)
				R.flag = jacobianECADD(&R, P, p, &R);
		}
	}
	jac2Aff(&R, p, kP);
}
void generalMontgomery(IN word* k, IN point* P, IN word* p, OUT point* kP)
//jacobian으로 하면 되지 않는다.
{
	int i, j;
	point R0, R1;
	R0.flag = 1;
	R1 = *P;

	for (i = SIZE - 1; i >= 0; i--)
	{
		for (j = 31; j >= 0; j--)
		{
			if ((k[i] >> j) & 0x1)
			{
				R0.flag = ECADD(&R0, &R1, p, &R0);
				R1.flag = ECDBL(&R1, p, &R1);
			}
			else
			{
				R1.flag = ECADD(&R0, &R1, p, &R1);
				R0.flag = ECDBL(&R0, p, &R0);
			}
		}
	}
	memcpy(kP, &R0, sizeof(point));
}
//word형 배열 k를 w bit 단위로 나누어서 k2에 저장한다.
//k2는 크기 ceil(256/w)의 배열
void windows(IN word* k, IN int w, IN word* p, OUT word* k2)
{
	int i, d, cnt;
	d = ceil(256 / w);
	word temp[SIZE];
	memcpy(temp, k, sizeof(word) * SIZE);
	
	for (i = 0; i < d; i++)
	{
		cnt = 0;
		k2[i] = temp[0] & ((1 << w) - 1);
		while (cnt++ != w)
			div2(temp);
	}
}
void precomputationWindow(IN point* P, IN int w, IN word* p, OUT point* pP)
{
	int i, j, num;
	num = (((long long)1 << w) - 1);
	word** k;
	k = (word**)malloc(sizeof(word*) * num);	//2^w-1개, 1*P~(2^w-1)P
	for (i = 0; i < num; i++)
		k[i] = (word*)malloc(sizeof(word) * SIZE);

	for (i = 0; i < num; i++)
		memset(k[i], 0, sizeof(word) * SIZE);
	
	k[0][0] = 1;
	for (i = 1; i < num; i++)
	{
		k[i][0] = k[i - 1][0] + 1;
	}

	for (i = 0; i < num; i++)
	{
		jacobianLtR(k[i], P, p, &pP[i]);
		//LtR(k[i], P, p, &pP[i]);
		pP[i].flag = 0;
	}
}
void fixedWindow(IN word* k2, IN point* pP, IN int w, IN word* p, OUT point* kP)
{
	int i, j, d;
	point R;
	R.flag = 1;
	d = ceil(256 / w);
	
	for (i = d - 1; i >= 0; i--)
	{
		for (j = 0; j < w; j++)
			R.flag = jacobianECDBL(&R, p, &R);
		if (k2[i])
			R.flag = jacobianECADD(&R, &pP[k2[i] - 1], p, &R);
	}
	jac2Aff(&R, p, kP);
}
int main()
{
	word p[SIZE] = { 0, }, k[SIZE] = { 0, };
	point P, Q;
	byte buf[100];

	FILE* fp1, * fp2, * fp3, * fp4;
	int i, j, num, d, w;
	ull cycles1 = 0, cycles2 = 0, cycles3 = 0, cycles4 = 0, cycles5 = 0, cycles6 = 0;

	p[7] = 0xffffffff; p[6] = 0x00000001; p[5] = 0x00000000; p[4] = 0x00000000;
	p[3] = 0x00000000; p[2] = 0xffffffff; p[1] = 0xffffffff; p[0] = 0xffffffff;

	P.x[7] = 0x6b17d1f2; P.x[6] = 0xe12c4247; P.x[5] = 0xf8bce6e5; P.x[4] = 0x63a440f2;
	P.x[3] = 0x77037d81; P.x[2] = 0x2deb33a0; P.x[1] = 0xf4a13945; P.x[0] = 0xd898c296;

	P.y[7] = 0x4fe342e2; P.y[6] = 0xfe1a7f9b; P.y[5] = 0x8ee7eb4a; P.y[4] = 0x7c0f9e16;
	P.y[3] = 0x2bce3357; P.y[2] = 0x6b315ece; P.y[1] = 0xcbb64068; P.y[0] = 0x37bf51f5;

	P.flag = 0;

	//Double_and_add, GeneralMontgomery, fixed_window
	fp1 = fopen("TV_Scalar.txt", "r");
	//fp2 = fopen("ret_DaA.txt", "w");
	//fp2 = fopen("ret_GM.txt", "w");
	//fp2 = fopen("ret_FW.txt", "w");

	fp2 = fopen("for_seminar.txt", "w");

	fp3 = fopen("test.txt", "r");
	fp4 = fopen("ret_test.txt", "w");

	/*printf("w값을 입력하시오 : ");
	scanf("%d", &w);*/
	w = 16;
	d = ceil(256 / w);
	num = (((long long)1 << w) - 2);
	word* k2;
	point* pP;
	k2 = (word*)malloc(sizeof(word) * d);
	pP = (word*)malloc(sizeof(point) * num);
	precomputationWindow(&P, w, p, pP);

	for (i = 0; i < 10000; i++)
	//for (i = 0; i < 3; i++)
	{
		//readData(k, fp3);
		readData(k, fp1);

		windows(k, w, p, k2);

		cycles1 = cpucycles();
		//Double_and_Add(k, &P, p, &Q);
		//generalMontgomery(k, &P, p, &Q);
		fixedWindow(k2, pP, w, p, &Q);
		cycles2 = cpucycles();

		cycles3 += (cycles2 - cycles1);

		for (j = SIZE - 1; j >= 0; j--)
		{
			fprintf(fp2, "%08X", Q.x[j]);
			//fprintf(fp4, "%08X", Q.x[j]);
		}
		fprintf(fp2, "\n");
		//fprintf(fp4, "\n");
		for (j = SIZE - 1; j >= 0; j--)
		{
			fprintf(fp2, "%08X", Q.y[j]);
			//fprintf(fp4, "%08X", Q.y[j]);
		}
		fprintf(fp2, "\n\n");
		//fprintf(fp4, "\n\n");

		fgets(buf, sizeof(buf), fp1);
		//fgets(buf, sizeof(buf), fp3);
	}
	printf("%10lld\n", cycles3 / 10000);

	return 0;
}
