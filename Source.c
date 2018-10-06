#include <stdio.h>
#include <math.h>
#include <stdlib.h>
// 7/5, 7.9, 3040000000, 0, 7/5, 11.37, 50000, 1179280000 
//5/3, 0.0002219, 3781200, -158700, 5/3, 0.00271, 10, 500000
// 5/3, 7.9, 3040000000, 0, 7/5, 11.37, 50000, 1179280000
/*long double mem = 65000;
const int memarr = 65000;*/

long double absl(long double r) {  //calculating absolute value of long double
	if (r < 0)
		r = r * (-1);
	return r;
}

long double dst(long double a, long double b) {
	long double y;
	y = (absl(a) + absl(b)) / 2;
	return y;
}

long double power(long double r, int n) { // calculating of long double * n
	if (n == 0) {
		r = 1;
		return r;
	}
	else {
		long double m;
		m = 1;
		for (int k = 1; k <= n; k++) {
			m = r * m;
		}
		return m;
	}
}

int comp(long double a, long double b) { // comparison of two long doubles
	if (a * b <= 0)
		return 1;
	else
		return 0;
}

//7/5, 7.9, 3040000000, 0, 7/5, 11.37, 50000, 1179280000 given values (1.5)


int main()
{

	long double gamma_0, rho_0, P_0, U_0, gamma_3, rho_3, U_3, P_3 = 0;
	long double g_0_1, g_0_2, g_1_1, g_1_2 = 0, A = 0, B = 0;
	printf("Enter the following data (divided by comma): gamma_0, rho_0, P_0, U_0, gamma_3, rho_3, U_3, P_3 \n");
	scanf_s("%Lf/%Lf, %Lf, %Lf, %Lf, %Lf/%Lf, %Lf, %Lf, %Lf", &g_0_1, &g_0_2, &rho_0, &P_0, &U_0, &g_1_1, &g_1_2, &rho_3, &U_3, &P_3);


	printf("%Lf %Lf %Lf %Lf\n", g_0_1, g_0_2, g_1_1, g_1_2);

	gamma_0 = (g_0_1 / g_0_2);
	gamma_3 = (g_1_1 / g_1_2);

	printf("%Lf %Lf %lf %lf %lf %lf %lf %lf\n", gamma_0, rho_0, P_0, U_0, gamma_3, rho_3, U_3, P_3);

	long double a_0, a_1, a_2, a_3, a_4, a_5, a_6 = 0;
	long double alpha_0, alpha_3 = 0;
	long double e_0, e_3, C_0, C_3, X = 0;
	X = P_3 / P_0;
	alpha_0 = (gamma_0 + 1) / (gamma_0 - 1);
	alpha_3 = (gamma_3 + 1) / (gamma_3 - 1);
	C_0 = sqrt(gamma_0 * P_0 / rho_0);
	C_3 = sqrt(gamma_3 * P_3 / rho_3);
	e_0 = (2 * C_0 * C_0) / (gamma_0 * (gamma_0 - 1) * (U_3 - U_0) * (U_3 - U_0));
	e_3 = (2 * C_3 * C_3) / (gamma_3 * (gamma_3 - 1) * (U_3 - U_0) * (U_3 - U_0));
	printf("%Lf %Lf %Lf %Lf %Lf %Lf\n", alpha_0, alpha_3, C_0, C_3, e_0, e_3);
	//calculating the coefficients

	a_0 = (alpha_0 * e_3 - alpha_3 * X * e_0) * (alpha_0 * e_3 - alpha_3 * X * e_0);
	a_1 = 2 * ((alpha_0 * e_3 - alpha_3 * X * e_0) * (e_3 * (1 - 2 * alpha_0 * X) - e_0 * X * (X - 2 * alpha_3))
		- alpha_0 * alpha_3 * X * (alpha_0 * e_3 + alpha_3 * X * e_0));
	a_2 = e_3 * e_3 * (6 * alpha_0 * alpha_0 * X * X - 8 * alpha_0 * X + 1) -
		2 * e_0 * e_3 * X * (alpha_0 * alpha_3 * (X * X + 4 * X + 1) - 2 * (X + 1) * (alpha_3 + alpha_0 * X) + X) +
		e_0 * e_0 * X * X * (6 * alpha_3 * alpha_3 - 8 * alpha_3 * X + X * X) + alpha_0 * alpha_0 * alpha_3 * alpha_3 * X * X -
		2 * alpha_0 * X * e_3 * (alpha_0 * X - 2 * alpha_0 * alpha_3 * X + 2 * alpha_3) -
		2 * alpha_3 * X * X * e_0 * (alpha_3 + 2 * alpha_0 * X - 2 * alpha_0 * alpha_3);
	a_3 = -2 * X * (2 * e_3 * e_3 * (alpha_0 * alpha_0 * X * X - 3 * alpha_0 * X + 1) +
		e_0 * e_3 * ((alpha_3 + alpha_0 * X) * (X * X + 4 * X + 1) - 2 * alpha_0 * alpha_3 * X * (X + 1) - 2 * X * (X + 1)) +
		2 * e_0 * e_0 * X * (X * X - 3 * alpha_3 * X + alpha_3 * alpha_3) - alpha_0 * alpha_3 * X * (alpha_0 * X + alpha_3) +
		e_3 * (alpha_0 * alpha_0 * alpha_3 * X * X - 2 * X * (2 * alpha_0 * alpha_3 + alpha_0 * alpha_0 * X) + (2 * alpha_0 * X + alpha_3)) +
		e_0 * X * (alpha_0 * alpha_3 * alpha_3 - 2 * alpha_3 * (alpha_3 + 2 * alpha_0 * X) + 2 * alpha_3 * X + alpha_0 * X * X));

	a_4 = X * X * (e_3 * e_3 * (alpha_0 * alpha_0 * X * X - 8 * alpha_0 * X + 6) -
		2 * e_0 * e_3 * (alpha_0 * alpha_3 * X - 2 * (X + 1) * (alpha_3 + alpha_0 * X) + X * X + 4 * X + 1) +
		e_0 * e_0 * (alpha_3 * alpha_3 - 8 * alpha_3 * X + 6 * X * X) + (alpha_3 * alpha_3 + 4 * alpha_0 * alpha_3 * X + alpha_0 * alpha_0 * X * X) -
		2 * e_3 * ((alpha_0 * alpha_0 * X + 2 * alpha_0 * alpha_3) * X - 2 * (2 * alpha_0 * X + alpha_3) + 1) -
		2 * e_0 * (alpha_3 * (2 * alpha_0 * X + alpha_3) - 2 * X * (2 * alpha_3 + alpha_0 * X) + X * X));
	a_5 = 2 * X * X * X * (e_3 * e_3 * (alpha_0 * X - 2) - e_0 * e_3 * (alpha_0 * X - 2 + alpha_3 - 2 * X) +
		e_0 * e_0 * (alpha_3 - 2 * X) + (alpha_3 + alpha_0 * X) -
		e_3 * (2 * alpha_0 * X + alpha_3 - 2) - e_0 * (2 * alpha_3 + alpha_0 * X - 2 * X));
	a_6 = X * X * X * X * ((e_3 - e_0) * (e_3 - e_0) + 1 - 2 * (e_3 + e_0));

	printf("Polynome coefficients: \n");
	printf("a_0 = %Lf\n", a_0);
	printf("a_1 = %Lf\n", a_1);
	printf("a_2 = %Lf\n", a_2);
	printf("a_3 = %Lf\n", a_3);
	printf("a_4 = %Lf\n", a_4);
	printf("a_5 = %Lf\n", a_5);
	printf("a_6 = %Lf\n", a_6);


	long double a[7]; // array for calculating A
	a[0] = a_0;
	a[1] = a_1;
	a[2] = a_2;
	a[3] = a_3;
	a[4] = a_4;
	a[5] = a_5;
	a[6] = a_6;

	long double b[7]; // array for calculating B
	b[0] = a_0;
	b[1] = a_1;
	b[2] = a_2;
	b[3] = a_3;
	b[4] = a_4;
	b[5] = a_5;
	b[6] = a_6;

	for (int q = 0; q < 7; q++) { // deriving absolute values of all coefficients
		a[q] = absl(a[q]);
		b[q] = absl(b[q]);
	}


	int i = 0, j = 0; // sorting to find A
	long double c = 0;
	for (i = 1; i < 6; i++) {
		for (j = 5; j >= i; j--) {
			if (a[j] > a[j + 1])
			{
				c = a[j];
				a[j] = a[j + 1];
				a[j + 1] = c;
			}
		}
	}

	for (i = 0; i < 5; i++) { // sorting to find B
		for (j = 4; j >= i; j--) {
			if (b[j] > b[j + 1])
			{
				c = b[j];
				b[j] = b[j + 1];
				b[j + 1] = c;
			}
		}
	}
	A = a[6];
	B = b[5];


	long double z1 = (absl(a_6)) / (B + absl(a_6)), z2 = 1 + A / absl(a_0); // finding the complex ring where the roots belong
																			//	printf("%Lf %Lf \n", z1, z2);

	long double seg = 0;
	const int segint = 60000;
	seg = (z2 - z1) / segint; // segment length
	long double point[60000][2];
	for (int m = 0; m <segint; m++) { // filling the array with coordinates of segments
		for (int n = 0; n < 2; n++) {
			if (n == 0) {
				point[m][n] = z1;
			}
			else {
				z1 = z1 + seg;
				point[m][n] = z1;
			}
		}
	}

	long double d[7]; // coefficients for deriving the polynome
	d[6] = a_0;
	d[5] = a_1;
	d[4] = a_2;
	d[3] = a_3;
	d[2] = a_4;
	d[1] = a_5;
	d[0] = a_6;

	long double answers[6][2]; // answers array
	for (int an1 = 0; an1 < 6; an1++) {
		for (int an2 = 0; an2 < 2; an2++) {
			answers[an1][an2] = 0;
		}
	}
	int r = 0;
	long double polinomeleft = 0;
	long double polinomeright = 0;
	int True = 0;
	for (int m = 0; m < 60000; m++) {
		for (int n = 0; n < 2; n++) {
			if (n == 0) {
				for (int p = 6; p >= 0; p--) {
					polinomeleft = d[p] * power(point[m][0], p) + polinomeleft; // definition of polynome on the left border of segment
				}
			}
			else {
				for (int wr = 6; wr >= 0; wr--) {
					polinomeright = d[wr] * power(point[m][1], wr) + polinomeright;// definition of polynome on the right border of segment
				}
			}

		}
		True = comp(polinomeleft, polinomeright);
		if (True == 1) {
			answers[r][0] = point[m][0]; // filling the answers array with the coordinates of segments where the roots belong
			answers[r][1] = point[m][1];
			r++;
			//					printf("wow\n");
			True = 0;
		}

		polinomeright = 0;
		polinomeleft = 0;
	}




	long double mo = -1;  // similar for the negatives
	z1 = z1*mo;
	z2 = z2*mo;
	for (int m = 0; m < 60000; m++) {
		for (int n = 0; n < 2; n++) {
			if (n == 0) {
				point[m][n] = z2;
			}
			else {
				z2 = z2 + seg;
				point[m][n] = z2;
			}
		}
	}

	for (int m = 0; m < 60000; m++) {
		for (int n = 0; n < 2; n++) {
			if (n == 0) {
				for (int p = 6; p >= 0; p--) {
					polinomeleft = d[p] * power(point[m][0], p) + polinomeleft;
				}
			}
			else {
				for (int wr = 6; wr >= 0; wr--) {
					polinomeright = d[wr] * power(point[m][1], wr) + polinomeright;
				}
			}

		}
		True = comp(polinomeleft, polinomeright);
		if ((True == 1) && (polinomeleft > 0) && (polinomeright > 0)) {  // the roots are positive
			answers[r][0] = point[m][0];
			answers[r][1] = point[m][1];
			r++;
			//			printf("wow\n");
			True = 0;
		}

		polinomeright = 0;
		polinomeleft = 0;
	}

	int TRUE = 0;
	printf("Localization intervals: \n");

	for (int r = 0; r < 6; r++) {
		for (int s = 0; s < 2; s++) {
			if (answers[r][s] != 0)
				TRUE = 1;
		}
		if (TRUE == 1)
			printf("[%Lf; %Lf]", answers[r][0], answers[r][1]); // printing out the answers array
		printf("\n");
		TRUE = 0;
	}

	long double e_left = 0;
	long double e_right = 0;
	long double middle = 0;
	long double ans[6];
	for (int yyy = 0; yyy < 6; yyy++) {
		ans[yyy] = 0;
	}
	long double polinomemiddle = 0;
	polinomeleft = 0;
	polinomeright = 0;

	printf("Clarified roots:\n");

	int rootn = 0;
	int increment = 6;
	long double roots[6];

	for (int y = 0; y < 6; y++) {
		if (answers[y][0] == 0 && answers[y][1] == 0)
			increment--;
	}


	for (i = 0; i < increment; i++) {
		e_left = answers[i][0];
		e_right = answers[i][1];
		middle = dst(e_left, e_right);
		for (int t = 0; t < 10000; t++) {
			middle = dst(e_left, e_right);
			for (int p = 6; p >= 0; p--) {
				polinomeleft = d[p] * power(e_left, p) + polinomeleft;
			}

			for (int p = 6; p >= 0; p--) {
				polinomeright = d[p] * power(e_right, p) + polinomeright;
			}
			for (int p = 6; p >= 0; p--) {
				polinomemiddle = d[p] * power(middle, p) + polinomemiddle;
			}
			if (polinomeleft == 0) {
				ans[i] = polinomeleft;
			}
			else if (polinomeright == 0) {
				ans[i] = polinomeright;
			}
			else if (polinomemiddle == 0) {
				ans[i] = polinomemiddle;
			}
			if (comp(polinomeleft, polinomemiddle) == 1) {
				e_right = middle;
			}
			else if (comp(polinomemiddle, polinomeright) == 1) {
				e_left = middle;
			}
			polinomemiddle = 0;
			polinomeleft = 0;
			polinomeright = 0;
		}
		if (ans[i] == 0) {
			printf("[%.10Lf; %.10Lf]\n", e_left, e_right);
		}
		else {
			printf("%Lf\n", ans[i]);
		}
		roots[i] = e_left;
		printf("%.10Lf\n", e_left);
	}

	printf("\n");

	long double D_0 = 0, D_3 = 0, rho_1 = 0, rho_2 = 0, U_1 = 0, P_1 = 0, P_2 = 0, Y = 0, U_2 = 0, a_cr, l_1, l_2;
	for (i = 0; i < increment; i++) {
		Y = roots[i];
		printf("Y = %.10Lf \n", Y);
		P_1 = Y * P_0;
		P_2 = P_1;
		rho_1 = rho_0 * ((gamma_0 - 1) + Y * (gamma_0 + 1)) / ((gamma_0 + 1) + Y * (gamma_0 - 1));
		U_1 = U_0 + sqrt((P_1 - P_0) * (rho_1 - rho_0) / (rho_1 * rho_0));
		rho_2 = rho_3 * ((gamma_3 + 1) * P_2 + (gamma_3 - 1) * P_3) / ((gamma_3 - 1) * P_2 + (gamma_3 + 1) * P_3);
		U_2 = U_3 - sqrt((P_2 - P_3) * (rho_2 - rho_3) / (rho_2 * rho_3));

		printf(" P_1 = %0.10Lf \n P_2 = %0.10Lf \n rho_1 = %0.10Lf \n rho_2 = %0.10Lf \n U_1 = %0.10Lf \n U_2 = %0.10Lf \n", P_1, P_2, rho_1, rho_2, U_1, U_2);

		D_0 = (rho_1 * U_1 - rho_0 * U_0) / (rho_1 - rho_0);

		a_cr = (P_1 - P_0) / (rho_1 - rho_0);
		l_1 = (U_0 - D_0) / sqrt(a_cr);
		l_2 = (U_1 - D_0) / sqrt(a_cr);
		printf("a_cr = %0.10Lf, l_1 = %0.10Lf, l_2 = %0.10Lf, l_1xl_2 = %0.10Lf \n", a_cr, l_1, l_2, l_1 * l_2);

		D_3 = (rho_2 * U_2 - rho_3 * U_3) / (rho_2 - rho_3);

		printf(" D_0 = %.10Lf \n D_3 = %.10Lf \n", D_0, D_3);

		printf("\n");
	}



	system("PAUSE");
	return 0;
}