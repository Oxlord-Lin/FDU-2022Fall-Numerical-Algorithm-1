#include <stdio.h>

int main() {
	long long n = 1 ;
	float sum1 = 0, sum2 = 0;
	while (1) {
		//single precision floating-point arithmetic
		sum1 = sum2;
		sum2 += (1 / (n + 0.0));
		if (sum1 == sum2)
			break;
		n++;
	}
	printf("\nÌø³öÑ­»·£¡\nsum=%.7lf,n=%lld", sum2, n - 1);
}