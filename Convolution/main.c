#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*
 * グラフ化した複素数をCSVファイルへ書き込む
 * ----
 * reh, imh : データバッファ（実部、虚部）
 * ref, imf : データバッファ（実部、虚部）
 * reg, img : データバッファ（実部、虚部）
 * n        : データ数
 */
int writeCsvFile(char *dst, double *real_f, double *real_g, double *real_h,
				 double *freq_f, double *freq_g, double *freq_h, int n)
{
	int ret = 0, i;
	FILE *fp = NULL;
	
	for (;;) {
		if ((fp = fopen(dst, "w")) == NULL) break;
		fprintf(fp, "f, g, f*g, F=FT(f), G=FT(g), h*g=iFT(FG)\n");
		for (i = 0; i < n; i++) {
			fprintf(fp, "%f, %f, %f, %f, %f, %f\n",
					real_f[i], real_g[i], real_h[i], freq_f[i], freq_g[i], freq_h[i]);
		}
		// 正常終了
		ret = 1;
		break;
	}
	if (fp) fclose(fp);
	return ret;
}

/*
 * 変数の値を入れ替える
 * ----
 * a, b : 入れ替える変数
 */
void swap(double *a, double *b)
{
	double tmp = *a;
	*a = *b;
	*b = tmp;
}

/*
 * 低周波と高周波のデータを入れ替える
 * ----
 * dat : データバッファ
 * n   : データ数
 */
void replace(double *dat, int n)
{
	int i, nh;
	
	nh = n / 2;
	for (i = 0; i < nh; i++)
		swap(&dat[i], &dat[i+nh]);
}

/*
 * 高速フーリエ変換を行う
 * ----
 * rev, imv : データバッファ（実部、虚部）
 * n        : データ数
 * inverse  : 逆変換フラグ（0:順変換、1:逆変換）
 * f()      : 成功フラグを返す（成功なら1、失敗なら0）
 */
int fft(double *rev, double *imv, int n, int inverse)
{
	int i, j, k, m, n1, n2, t;
	double theta, wr, wi, dr, di;
	
	// バタフライ演算
	theta = ((inverse)? M_PI: -M_PI) / (double)n;
	for (i = n, m = 0; 2 <= i; i >>= 1, m++)
		;
	n2 = n;
	for (k = 0; k < m; k++) {
		theta *= 2.0;
		n1 = n2;
		n2 = n1 / 2;
		for (j = 0; j < n2; j++) {
			wr = cos(theta*j);
			wi = sin(theta*j);
			for (i = j; i < n; i+=n1) {
				t = i + n2;
				dr = rev[i] - rev[t];
				di = imv[i] - imv[t];
				rev[i] += rev[t];
				imv[i] += imv[t];
				rev[t] = dr*wr - di*wi;
				imv[t] = di*wr + dr*wi;
			}
		}
	}
	// データを入れ替える
	for (i = j = 0; i < n-1; i++, j+=n2) {
		if (i < j) {
			swap(&rev[i], &rev[j]);
			swap(&imv[i], &imv[j]);
		}
		for (n2 = n/2; n2 <= j; j-=n2, n2>>=1)
			;
	}
	// 順変換の場合は1/nする
	if (! inverse) {
		for (i = 0; i < n; i++) {
			rev[i] /= (double)n;
			imv[i] /= (double)n;
		}
	}
	return 1;
}

/*
 * 値をリミットを制御する
 * ----
 * value    : 値
 * min, max : 最小値、最大値
 * f()      : 正規化された値
 */
int normalize(int value, int min, int max)
{
	if (value < min) value = min;
	else if (max < value) value = max;
	return value;
}

/*
 * Shepp-Logan フィルタを作成する
 * ----
 * rev, imv : データバッファ（実部、虚部）
 * n        : データ数
 * a        : サンプリング間隔
 */
void makeSheppLoganFilter(double *rev, double *imv, int n)
{
	int i, k;
	
	for (i = 0, k = n/2; i < n; i++, k--) {
		rev[i] = 2.0 / (M_PI*M_PI*(1.0-4.0*k*k));
		imv[i] = 0.0;
	}
}

/*
 * パルスデータを作成する
 * ----
 * rev, imv : データバッファ（実部、虚部）
 * n        : データ数
 * peak     : ピーク値
 */
void makePulse(double *rev, double *imv, int n, int peak)
{
	int i;
	
	for (i = 0; i < n; i++)
		rev[i] = imv[i] = 0.0;
	for (i = n/2-n/16; i < n/2+n/16; i++)
		rev[i] = (double)peak;
}

/*
 * 複素数の絶対値を計算する
 * ----
 * rev, imv : データバッファ（実部、虚部）
 * n        : データ数
 */
void zabs(double *rev, double *imv, int n)
{
	int i;
	
	for (i = 0; i < n; i++) {
		rev[i] = sqrt(rev[i]*rev[i]+imv[i]*imv[i]);
		imv[i] = 0.0;
	}
}

/*
 * 複素数のスカラー倍を計算する
 * ----
 * rev, imv : データバッファ（実部、虚部）
 * n        : データ数
 * k        : スカラー
 */
void zscalar(double *rev, double *imv, int n, int k)
{
	int i;
	
	for (i = 0; i < n; i++) {
		rev[i] *= (double)k;
		imv[i] *= (double)k;
	}
}

/*
 * 畳み込み積分を計算をする
 * ----
 * reh, imh : データバッファ（実部、虚部）
 * ref, imf : フィルタバッファ（実部、虚部）
 * reg, img : データバッファ（実部、虚部）
 * n        : データ数
 * m        : フィルタに使用するデータ数
 */
void zconvo(double *reh, double *imh,
			double *ref, double *imf, double *reg, double *img, int n, int m)
{
	int i, j, ofsf, ofsg, indexf, indexg;
	
	// 合計
	ofsf = n/2 - m/2;
	double sum = 0.0;
	for (j = 0; j < m; j++) {
		indexf = normalize(ofsf+j, 0, n-1);
		sum += ref[indexf];
	}
	printf("sum = %f\n", sum);
	// 畳み込み積分（コンボリューション）を計算する
	ofsf = n/2 - m/2;
	for (i = 0; i < n; i++) {
		reh[i] = imh[i] = 0.0;
		ofsg = i - m/2;
		for (j = 0; j < m; j++) {
			indexf = normalize(ofsf+j, 0, n-1);
			indexg = normalize(ofsg+j, 0, n-1);
			reh[i] += (double)(ref[indexf]*reg[indexg]);
		}
		reh[i] /= sum;
	}
	
}

/*
 * 複素数の乗算を計算をする
 * ----
 * reh, imh : データバッファ（実部、虚部）
 * ref, imf : データバッファ（実部、虚部）
 * reg, img : データバッファ（実部、虚部）
 * n        : データ数
 */
void zmul(double *reh, double *imh,
		  double *ref, double *imf, double *reg, double *img, int n)
{
	int i;
	
	for (i = 0; i < n; i++) {
		reh[i] = ref[i]*reg[i] - imf[i]*img[i];
		imh[i] = ref[i]*img[i] + imf[i]*reg[i];
	}
}

/*
 * 複素数データをクラフ化のために保存する
 * ----
 * dat        : データバッファ
 * rev, imv   : データバッファ（実部、虚部）
 * n          : データ数
 * freq_space : 領域の種類（0:実領域、1:周波数領域）
 */
void recResult(double *dat, double *rev, double *imv, int n, int freq_space)
{
	int i;
	
	if (freq_space) {
		for (i = 0; i < n; i++)
			dat[i] = sqrt(rev[i]*rev[i]+imv[i]*imv[i]);
		replace(dat, n);
	} else {
		for (i = 0; i < n; i++)
			dat[i] = rev[i];
	}
}

/*
 * 実領域と周波数領域でフィルタをかけて結果をCSVファイルに出力する
 * ----
 * dst : 出力ファイル名
 * f() : 成功フラグを返す（成功なら1、失敗なら0）
 */
int process(char *dst)
{
	enum {real_f, real_g, real_h, freq_f, freq_g, freq_h,
		ref, imf, reg, img, reh, imh, BUFN};
	enum {forward, inverse, real_space=0, freq_space};
	const int n = 64, m = 13;
	int ret = 1, i;
	double *buf[BUFN];
	
	for (i = 0; i < BUFN; i++) buf[i] = NULL;
	for (;;) {
		// メモリを確保する
		for (i = 0; i< BUFN; i++)
			if ((buf[i] = (double *)malloc(sizeof(double)*n)) == NULL) {ret=0; break;}
		if (ret == 0) break;
		// fを作成する
		makeSheppLoganFilter(buf[ref], buf[imf], n);
		recResult(buf[real_f], buf[ref], buf[imf], n, real_space);
		// gを作成する
		makePulse(buf[reg], buf[img], n, 30);
		recResult(buf[real_g], buf[reg], buf[img], n, real_space);
		// f*gを計算する
		zconvo(buf[reh], buf[imh], buf[ref], buf[imf], buf[reg], buf[img], n, m);
		recResult(buf[real_h], buf[reh], buf[imh], n, real_space);
		// Fを準備する
		fft(buf[ref], buf[imf], n, forward);
		zabs(buf[ref], buf[imf], n);
		recResult(buf[freq_f], buf[ref], buf[imf], n, freq_space);
		// Gを準備する
		fft(buf[reg], buf[img], n, forward);
		recResult(buf[freq_g], buf[reg], buf[img], n, freq_space);
		// FGを計算する
		zmul(buf[reh], buf[imh], buf[ref], buf[imf], buf[reg], buf[img], n);
		fft(buf[reh], buf[imh], n, inverse);
		zscalar(buf[reh], buf[imh], n, n);
		recResult(buf[freq_h], buf[reh], buf[imh], n, real_space);
		// ファイルを出力する
		ret = writeCsvFile(dst, buf[real_f], buf[real_g], buf[real_h],
						   buf[freq_f], buf[freq_g], buf[freq_h], n);
		break;
	}
	// メモリを解放する
	for (i = BUFN-1; 0 <= i; i--)
		if (buf[i]) free(buf[i]);
	return ret;
}

int main(void)
{
	char dst[256];
	
	printf("出力ファイル： ");	scanf("%s", dst);
	printf((process(dst))? ">成功しました。\n": ">失敗しました。\n");
	return 0;
}
