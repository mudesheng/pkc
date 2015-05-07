package discretecosinetransform

// #cgo CFLAGS: -I/usr/local/include -I/usr/include/lapacke/ -std=c99
// #cgo LDFLAGS: -L ./ -lm -lcv -lcxcore -ljpeg -lpng -lhighgui -llapack -llapacke
// #include <cv.h>
// #include <stdio.h>
// #include <math.h>
// #include <highgui.h>
// #include <cxcore.h>
// #include <lapacke.h>
/*


typedef struct {
  double *s;
  double *u;
  double *vt;
} SV ;

typedef struct {
  double *cdctarr;
  CvMat *dst;
} DCTRet;

void jpgcvDCT() {
  IplImage *input, *output, *b;
  input = cvLoadImage("dct_orig.png",CV_LOAD_IMAGE_COLOR);
  int width = input->width;
  int height = input->height;
  // b = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U,1);
  b = cvCreateImage(cvSize(width, height), CV_32FC1,1);
  IplImage *b_dct;
  b_dct=cvCreateImage(cvSize(width,height),CV_32FC1,1);
  cvDCT(b,b_dct,0);  // doubt??
  cvShowImage("blue",b_dct);
  cvSaveImage("/db/goworkspace/src/pkc/dct.jpg", b_dct);
}

double* NewdoubleArray(int length) {
  double *arr = calloc(length, sizeof(double));
  return arr;
}

inline void SetdoubleArray(double *arr, int pos, double val) {
  arr[pos] = val;
}

// Auxiliary routine: printing a matrix
void print_matrix( char* desc, int m, int n, double* a, int lda  ) {
	int i, j;
	printf( "\n %s\n", desc  );
	for( i = 0; i < m; i++  ) {
		for( j = 0; j < n; j++  ) printf( " %6.2f", a[i*lda+j]);
		printf( "\n"  );
	}
}

void print_cvMat(CvMat *mat) {
  printf("[print_cvMat]cvMat:\n");
  for(int i = 0; i <mat->rows;i++) {
    for(int j = 0; j < mat->cols; j++) {
      printf(" %6.2f", cvmGet(mat, i, j));
    }
    printf("\n");
  }
}

DCTRet ReadsDCT(double *data, const int length) {
  // printf("data pointer address: %p\t", data);
  // printf("length: %d, clang arr:\t", length);
  // printf("%d:%f, %d:%f\n", 0, data[0], length-1, data[length-1]);
  // for( int i = 0; i < length; i++) {
  //   printf("%f ");
  // }
  // printf("\n");
  double row = sqrt(length);
  if(row * row != length  && (int)row % 2 == 0) {
    fprintf(stderr, "[ReadsDCT]the input data matrix rows and columns must be equaled and even, exit...\n");
      exit(1);
  }

  CvMat src64 ;
  CvMat *dst64 = cvCreateMat(row, row, CV_64FC1);
  cvInitMatHeader(&src64, row,row ,CV_64FC1, data, CV_AUTOSTEP);
  cvDCT(&src64,dst64,CV_DXT_FORWARD);
  // print_cvMat(dst64);
  // print_matrix("dst64", dst64->rows, dst64->cols, (double*)dst64->data.db, dst64->cols);
  // free work
  cvReleaseData(&src64);
  // free(data);

  DCTRet dctRet;
  dctRet.cdctarr = (double*)dst64->data.db;
  dctRet.dst = dst64;
  return dctRet;

}

SV PSVD() {
  pdgesvd('V','N',m, n, dataMat, 1,1,)
}

SV SVD(double *dataMat, const int row, const int column) {
	// Locals
	int M = row, N = column;
	int LDA = N, LDU = M, LDVT = N;
	int m = M, n = N, lda = LDA, ldu = LDU, ldvt = LDVT,info;

	int min = row > column ? column :row;
	// double superb[min-1];
	double *superb = calloc(min-1, sizeof(double));
	// double s[N], u[LDU *M], vt[LDVT *N];
	double *s = calloc(N, sizeof(double));
	// double *u = calloc(LDU * M, sizeof(double));
	double *u = calloc(LDU * M, sizeof(double));
	double *vt = calloc(LDVT *N, sizeof(double));

	// Executable statements
	printf( "[C.SVD] DGESVD Example Program Results\n");
	fflush(stdout);

	// Compute SVD
	// info = LAPACKE_dgesvd( LAPACK_ROW_MAJOR,'A', 'A', m, n, dataMat, lda, s, u, ldu, vt, ldvt,superb);
	info = LAPACKE_dgesvd( LAPACK_ROW_MAJOR,'O', 'N', m, n, dataMat, lda, s, u, ldu, NULL, ldvt,superb);
	// dgesvd("All", "All", &m, &n, dataMat, &lda, s, u, &ldu, vt, &ldvt, wkopt, &lwork);
	// lwork = (int)wkopt;
	// work = (double*)malloc(lwork * sizeof(double));
	// LAPACKE_sgesvd("All", "All", &m, &n, dataMat, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);
	// // check for convergence
	if(info > 0) {
	  fprintf(stderr, "[__FUNC__] The algorithm computing SVD failed to converge.\n");
	  fflush(stderr);
	  exit(1);
	}
	// Print singular values
	print_matrix( "Singular values", 1, n, s, 1);
	// Print left singular vectors
	print_matrix( "Left singular vectors (stored columnwise)", m, n, u, ldu  );
	// Print right singular vectors
	print_matrix( "Right singular vectors (stored rowwise)", n, n, vt, ldvt  );
	// Free workspace
	free(superb);
	// free(s); free(u); free(vt);
	SV sv;
	sv.s = s; sv.u = u; sv.vt = vt;

	return sv;
}

  // IplImage *image=0;
  // IplImage *image1=0;
  //  IplImage *image2=0;
  // IplImage *image3=0;
  // CvMat* src32 = NULL;
  // CvMat arr;
  //       CvMat* dst=NULL;
  //       CvMat* idst=NULL;
	// double data[] = {1,2,3,4,\
			    // 6,7,8,10,\
			    // 11,12,14,15,\
			    // 16,17,19,20,\
			    // 21,22,24,25};

	// image=cvLoadImage("timeline.jpg",CV_LOAD_IMAGE_COLOR);

	  // image1=cvCreateImage(cvSize(image->width,image->height),IPL_DEPTH_32F,1);
	    // image2=cvCreateImage(cvSize(image->width, image->height),IPL_DEPTH_32F,1);
	      // image3=cvCreateImage(cvSize(image->width, image->height),IPL_DEPTH_32F,1);
	      // src32=cvCreateMat(image->height, image->width,CV_32FC1);
	      // dst=cvCreateMat(image->height, image->width, CV_64FC1);
	      // idst=cvCreateMat(image->height, image->width, CV_64FC1);

	      // // cvCvtColor(image,image1,CV_YCrCb2BGR);
	      // cvInitMatHeader(&arr, 4,4,CV_64FC1, data, 0);
	      // // cvCvtColor(image,image1,CV_RGB2Luv);
	         // // cvConvertScale(image1, src32,1,0);
		    // cvDCT(&arr,dst,CV_DXT_FORWARD);
		    // // cvDCT(im8,image2,CV_DXT_FORWARD);
		        // // cvConvertImage(image2,dst,0);
			// printf("dst:\n");
			// int i;
			// int h =4, w = 4;
			// for(i=0;i<h;i++)
			// {
			  // int j;
			   // for(j=0;j<w;j++)
			   // {
			       // printf("%5d",abs(cvmGet(dst,i,j)));

			   // }
			    // printf("\n");

			// }
			// cvDCT(dst,idst,CV_DXT_INVERSE);
			// printf("idst:\n");
			// for(i=0;i<h;i++)
			// {
			  // int j;
			   // for(j=0;j<w;j++)
			   // {
			       // printf("%5d",abs(cvmGet(idst,i,j)));

			   // }
			    // printf("\n");

			// }
			// // cvNamedWindow("jpeg",1);
			// //    cvNamedWindow("gray",1);
			// //       cvNamedWindow("dct",1);
			// //       cvNamedWindow("idct",1);
			          // // cvShowImage("jpeg",image);
				     // // cvShowImage("gray",image1);
				        // // cvShowImage("dct",image2);
					   // // cvShowImage("idct",image3);
					   // // cvSaveImage("dct.jpg", idst);
					       // // cvWaitKey(0);
					           // cvReleaseImage(&image);
						   // cvReleaseImage(&image1);
						   // cvReleaseImage(&image2);
						   // cvReleaseImage(&image3);
						   // cvReleaseMat(&src32);
						   // cvReleaseMat(&dst);

// }
*/
import "C"

import (
	"compress/gzip"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"sort"
	"unsafe"

	"code.google.com/p/biogo/alphabet"
	"code.google.com/p/biogo/io/seqio/fasta"
	"code.google.com/p/biogo/seq/linear"
	"github.com/jwaldrip/odin/cli"
	smatrix "github.com/skelterjohn/go.matrix"
)

type DCTArr struct {
	dct []float64
}

type SvdVect struct {
	ID     string
	LeftSV []float64
}

var Deviation float64 = 0.05

type SvdVectArr []SvdVect

func (sva SvdVectArr) Len() int {
	return len(sva)
}

func (sva SvdVectArr) Less(i, j int) bool {
	for x, y := 0, 0; x < len(sva[i].LeftSV) && y < len(sva[j].LeftSV); x, y = x+1, y+1 {
		if math.Abs(sva[i].LeftSV[x]-sva[j].LeftSV[y]) < sva[i].LeftSV[x]*Deviation {
			continue
		} else {
			return sva[i].LeftSV[x] < sva[j].LeftSV[y]
		}
	}

	return len(sva[i].LeftSV) < len(sva[j].LeftSV)
}

func (sva SvdVectArr) Swap(i, j int) {
	sva[i], sva[j] = sva[j], sva[i]
}

var DNA2Num = [5]C.double{0, 63, 127, 191, 255}

func DNA2Pixel(tl alphabet.Letters) (srcarr []C.double) {
	for _, e := range tl {
		switch e {
		case 'A':
			srcarr = append(srcarr, DNA2Num[1])
		case 'C':
			srcarr = append(srcarr, DNA2Num[2])
		case 'G':
			srcarr = append(srcarr, DNA2Num[3])
		case 'T':
			srcarr = append(srcarr, DNA2Num[4])
		default:
			srcarr = append(srcarr, DNA2Num[0])
		}
	}

	return
}

func AdjustSrcArr(dctarr []C.double) []C.double {
	sq := math.Sqrt(float64(len(dctarr)))
	if sq-math.Floor(sq) != 0 || int(sq)%2 != 0 {
		size := int(math.Floor(sq)) + 1
		if size%2 != 0 {
			size += 1
		}
		size *= size
		for i := len(dctarr); i < size; i++ {
			dctarr = append(dctarr, DNA2Num[0])
		}
	}

	return dctarr
}

func Gofloat64Arr2CdoubleArr(srcarr []float64, csrcarr unsafe.Pointer) {
	p := uintptr(csrcarr)
	for _, e := range srcarr {
		// *(*C.double)(unsafe.Pointer(&p)) = C.double(e)
		*(*C.double)(unsafe.Pointer(&p)) = C.double(e)
		fmt.Printf("%v\t", p)
		fmt.Printf("%f\t", *(*float64)(csrcarr))
		p += unsafe.Sizeof(e)
		fmt.Printf("%v\t", p)
	}
	fmt.Printf("%f\t", *(*float64)(csrcarr))
}

func CdoubleArr2Gofloat64Arr(cdctarr unsafe.Pointer, size int) (dctarr []float64) {
	p := uintptr(cdctarr)

	for i := 0; i < size; i++ {
		j := *(*float64)(unsafe.Pointer(p))
		dctarr = append(dctarr, j)
		p += unsafe.Sizeof(j)
	}

	return dctarr
}

func PrintCdctarr(cdctarr *C.double, length int) {
	p := uintptr(unsafe.Pointer(cdctarr))
	row := int(math.Sqrt(float64(length)))
	for i := 0; i < row; i++ {
		for j := 0; j < row; j++ {
			x := *(*float64)(unsafe.Pointer(p))
			fmt.Printf("%f ", x)
			p += unsafe.Sizeof(x)
		}
		fmt.Println()
	}
}

// get the most high power(left-up) square matrix
func GetPrimeDCTArr(cdctarr *C.double, length int, rows, cols int) (transdct []C.double) {
	p := uintptr(unsafe.Pointer(cdctarr))
	if rows*cols >= length {
		log.Fatalf("[GetPrimeDCTArr] Get size: %dX%d matrix bigger than length:%d\n", rows, cols, length)
	}

	x := *(*C.double)(unsafe.Pointer(p))
	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			tp := p + uintptr(i*cols+j)*unsafe.Sizeof(x)
			x := *(*C.double)(unsafe.Pointer(tp))
			transdct = append(transdct, x)
		}
	}

	return transdct
}

func DCT(cmd cli.Command) {
	// C.jpgcvDCT()
	// C.ReadsDCT()
	// double* ReadsDCT(double data[], const int length) {
	faName := cmd.Parent().Flag("f").String()
	fafp, err := os.Open(faName)
	if err != nil {
		log.Fatal("[DCT]open fa input file failed")
	}
	defer fafp.Close()
	fagzip, err := gzip.NewReader(fafp)
	far := fasta.NewReader(fagzip, linear.NewSeq("", nil, alphabet.DNA))
	var dctMatArr []C.double
	var rows, cols int = 30, 30
	var svdVectArr SvdVectArr

	for {
		if s, err := far.Read(); err != nil {
			if err == io.EOF {
				break
			} else {
				log.Fatalf("failed read %v:%s\n", faName, err)
			}
		} else {
			// fmt.Printf("ID: %v\n", s.CloneAnnotation().ID)
			t := s.(*linear.Seq)
			tl := t.Slice().(alphabet.Letters)
			// skip the length smaller than rows*cols
			if len(tl) < rows*cols {
				continue
			}
			var svdVect SvdVect
			svdVect.ID = s.CloneAnnotation().ID
			svdVectArr = append(svdVectArr, svdVect)
			srcarr := DNA2Pixel(tl)
			srcarr = AdjustSrcArr(srcarr)
			fmt.Printf("array len: %d, srcarr: %d:%f, %d:%f\n", len(srcarr), 0, srcarr[0], len(srcarr)-1, srcarr[len(srcarr)-1])
			// csrcarr := C.NewdoubleArray(C.int(len(srcarr)))
			// defer C.free(unsafe.Pointer(csrcarr))
			// fmt.Printf("Transform before:%v\n", csrcarr)
			// Gofloat64Arr2CdoubleArr(srcarr, unsafe.Pointer(csrcarr))
			// fmt.Printf("Transform after:%v\n", csrcarr)
			dctRet := C.ReadsDCT(&srcarr[0], C.int(len(srcarr)))
			cdctarr := dctRet.cdctarr
			dst := dctRet.dst
			// defer C.free(unsafe.Pointer(cdctarr))
			// dctarr := CdoubleArr2Gofloat64Arr(unsafe.Pointer(cdctarr), len(srcarr))
			// fmt.Printf("dctarr: %v\n", cdctarr)
			// PrintCdctarr(cdctarr, len(srcarr))
			transdct := GetPrimeDCTArr(cdctarr, len(srcarr), rows, cols)
			if len(transdct) != rows*cols {
				log.Fatalf("len(transdct):%d	!= rows*cols:%d\n", len(transdct), rows*cols)
			}
			dctMatArr = append(dctMatArr, transdct...)
			// var a DCTArr
			// a.dct = dctarr
			// dctMat = append(dctMat, a)
			// C.free(unsafe.Pointer(cdctarr))
			C.cvReleaseMat(&dst)

		}
	}

	// compute SVD matrix

	vectNum := len(dctMatArr) / (rows * cols)
	fmt.Printf("[SVD]len(dctMatArr): %d, vectNum: %v, len(svdVectArr): %d\n", len(dctMatArr), vectNum, len(svdVectArr))
	svM := C.SVD(&dctMatArr[0], C.int(vectNum), C.int(rows*cols))
	C.free(unsafe.Pointer(svM.vt))
	fmt.Printf("\n[SVD]svM: %v, svM.s: %v, svM.u: %v\n", svM, svM.s, svM.u)
	os.Stdout.Sync()
	singluarV := CdoubleArr2Gofloat64Arr(unsafe.Pointer(svM.s), rows*cols)
	C.free(unsafe.Pointer(svM.s))
	fmt.Printf("[SVD]singluarV: %v\n", singluarV)
	os.Stdout.Sync()
	// stored by columnwise
	e := *(*float64)(unsafe.Pointer(svM.u))
	for i := 0; i < vectNum; i++ {
		var leftVect []float64
		for j := 0; j < rows*cols; j++ {
			p := uintptr(unsafe.Pointer(svM.u)) + unsafe.Sizeof(e)*uintptr(i*vectNum+j)
			v := *(*float64)(unsafe.Pointer(p))
			leftVect = append(leftVect, v)
		}
		// leftVect := CdoubleArr2Gofloat64Arr(unsafe.Pointer(p), rows*cols)
		svdVectArr[i].LeftSV = leftVect
		// fmt.Printf("[SVD]LeftSV: %v\n", svdVectArr[i])
	}
	C.free(unsafe.Pointer(svM.u))

	// sort left Singluar Vectors
	fmt.Printf("After sort Left Singluar Vectors\n")
	os.Stdout.Sync()
	sort.Sort(svdVectArr)
	for _, item := range svdVectArr {
		fmt.Printf("%v: %v\n", item.ID, item.LeftSV)
	}
}
