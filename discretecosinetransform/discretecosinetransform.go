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

float* NewFloatArray(int length) {
  float *arr = calloc(length, sizeof(float));
  return arr;
}

inline void SetFloatArray(float *arr, int pos, float val) {
  arr[pos] = val;
}

float* ReadsDCT(float *data, const int length) {
  printf("data pointer address: %p\t", data);
  printf("length: %d, clang arr:\t", length);
  printf("%d:%f, %d:%f\n", 0, data[0], length-1, data[length-1]);
  // for( int i = 0; i < length; i++) {
  //   printf("%f ");
  // }
  printf("\n");
  float row = sqrt(length);
  if(row * row != length  && (int)row % 2 == 0) {
    fprintf(stderr, "[ReadsDCT]the input data matrix rows and columns must be equaled and even, exit...\n");
      exit(1);
  }

  CvMat src32 ;
  CvMat *dst32 = cvCreateMat(row, row, CV_32FC1);
  cvInitMatHeader(&src32, row,row ,CV_32FC1, data, CV_AUTOSTEP);
  cvDCT(&src32,dst32,CV_DXT_FORWARD);

  return (float*)dst32->data.fl;

}

// Auxiliary routine: printing a matrix
void print_matrix( char* desc, int m, int n, float* a, int lda  ) {
	int i, j;
	printf( "\n %s\n", desc  );
	for( i = 0; i < m; i++  ) {
		for( j = 0; j < n; j++  ) printf( " %6.2f", a[i*lda+j]);
		printf( "\n"  );
	}
}

void SVD(float *dataMat, const int row, const int column) {
	// Locals
	int M = row, N = column;
	int LDA = N, LDU = M, LDVT = N;
	int m = M, n = N, lda = LDA, ldu = LDU, ldvt = LDVT,info;

	int min = row > column ? column :row;
	float superb[min-1];
	float s[N], u[LDU *M], vt[LDVT *N];

	// Executable statements
	printf( " DGESVD Example Program Results\n"  );

	// Compute SVD
	info = LAPACKE_sgesvd( LAPACK_ROW_MAJOR,'A', 'A', m, n, dataMat, lda, s, u, ldu, vt, ldvt,superb);
	// dgesvd("All", "All", &m, &n, dataMat, &lda, s, u, &ldu, vt, &ldvt, wkopt, &lwork);
	// lwork = (int)wkopt;
	// work = (float*)malloc(lwork * sizeof(float));
	// LAPACKE_sgesvd("All", "All", &m, &n, dataMat, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);
	// // check for convergence
	if(info > 0) {
	  fprintf(stderr, "[__FUNC__] The algorithm computing SVD failed to converge.\n");
	  exit(1);
	}
	// Print singular values
	print_matrix( "Singular values", 1, n, s, 1  );
	// Print left singular vectors
	print_matrix( "Left singular vectors (stored columnwise)", m, n, u, ldu  );
	// Print right singular vectors
	print_matrix( "Right singular vectors (stored rowwise)", n, n, vt, ldvt  );
	// Free workspace
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
	"unsafe"

	"code.google.com/p/biogo/alphabet"
	"code.google.com/p/biogo/io/seqio/fasta"
	"code.google.com/p/biogo/seq/linear"
	"github.com/jwaldrip/odin/cli"
)

type DCTArr struct {
	dct []float64
}

var DNA2Num = [5]C.float{0, 63, 127, 191, 255}

func DNA2Pixel(tl alphabet.Letters) (srcarr []C.float) {
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

func AdjustSrcArr(dctarr []C.float) []C.float {
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

func GoFloat32Arr2CFloatArr(srcarr []float32, csrcarr unsafe.Pointer) {
	p := uintptr(csrcarr)
	for _, e := range srcarr {
		// *(*C.float)(unsafe.Pointer(&p)) = C.float(e)
		*(*C.float)(unsafe.Pointer(&p)) = C.float(e)
		fmt.Printf("%v\t", p)
		fmt.Printf("%f\t", *(*float32)(csrcarr))
		p += unsafe.Sizeof(e)
		fmt.Printf("%v\t", p)
	}
	fmt.Printf("%f\t", *(*float32)(csrcarr))
}

func CFloatArr2GoFloat32Arr(cdctarr unsafe.Pointer, size int) (dctarr []float32) {
	p := uintptr(cdctarr)

	for i := 0; i < size; i++ {
		j := *(*float32)(unsafe.Pointer(p))
		dctarr = append(dctarr, j)
		p += unsafe.Sizeof(j)
	}

	return dctarr
}

func PrintCdctarr(cdctarr *C.float, length int) {
	p := uintptr(unsafe.Pointer(cdctarr))
	row := int(math.Sqrt(float64(length)))
	for i := 0; i < row; i++ {
		for j := 0; j < row; j++ {
			x := *(*float32)(unsafe.Pointer(p))
			fmt.Printf("%f ", x)
			p += unsafe.Sizeof(x)
		}
		fmt.Println()
	}
}

func GetPrimeDCTArr(cdctarr *C.float, length int, retsize int) (transdct []C.float) {
	p := uintptr(unsafe.Pointer(cdctarr))
	if retsize >= length {
		log.Fatalf("[GetPrimeDCTArr] Get size:%d bigger than length:%d\n", retsize, length)
	}
	row := int(math.Sqrt(float64(length)))
	// i note row , j note columu
	x := *(*C.float)(unsafe.Pointer(p))
	for j := 0; j < row; j++ {
		for i, z := 0, j; i < row && z >= 0; i, z = i+1, z-1 {
			tp := p + uintptr(i*row+z)*unsafe.Sizeof(x)
			x := *(*C.float)(unsafe.Pointer(tp))
			transdct = append(transdct, x)
			if len(transdct) >= retsize {
				return transdct
			}
		}
	}

	return transdct
}

func DCT(cmd cli.Command) {
	// C.jpgcvDCT()
	// C.ReadsDCT()
	// float* ReadsDCT(float data[], const int length) {
	faName := cmd.Parent().Flag("f").String()
	fafp, err := os.Open(faName)
	if err != nil {
		log.Fatal("[DCT]open fa input file failed")
	}
	defer fafp.Close()
	fagzip, err := gzip.NewReader(fafp)
	far := fasta.NewReader(fagzip, linear.NewSeq("", nil, alphabet.DNA))
	var dctMatArr []C.float
	var inputDCTArrLen int = 1024

	for {
		if s, err := far.Read(); err != nil {
			if err == io.EOF {
				break
			} else {
				log.Fatalf("failed read %v:%s\n", faName, err)
			}
		} else {
			t := s.(*linear.Seq)
			tl := t.Slice().(alphabet.Letters)
			// skip the length smaller than inputDCTArrLen
			if len(tl) < inputDCTArrLen {
				continue
			}
			srcarr := DNA2Pixel(tl)
			srcarr = AdjustSrcArr(srcarr)
			fmt.Printf("array len: %d, srcarr: %d:%f, %d:%f\n", len(srcarr), 0, srcarr[0], len(srcarr)-1, srcarr[len(srcarr)-1])
			// csrcarr := C.NewFloatArray(C.int(len(srcarr)))
			// defer C.free(unsafe.Pointer(csrcarr))
			// fmt.Printf("Transform before:%v\n", csrcarr)
			// GoFloat32Arr2CFloatArr(srcarr, unsafe.Pointer(csrcarr))
			// fmt.Printf("Transform after:%v\n", csrcarr)
			cdctarr := C.ReadsDCT(&srcarr[0], C.int(len(srcarr)))
			// defer C.free(unsafe.Pointer(cdctarr))
			// dctarr := CFloatArr2GoFloat32Arr(unsafe.Pointer(cdctarr), len(srcarr))
			// fmt.Printf("dctarr: %v\n", cdctarr)
			PrintCdctarr(cdctarr, len(srcarr))
			transdct := GetPrimeDCTArr(cdctarr, len(srcarr), inputDCTArrLen)
			dctMatArr = append(dctMatArr, transdct...)
			// var a DCTArr
			// a.dct = dctarr
			// dctMat = append(dctMat, a)

		}
	}

	// compute SVD matrix
	row := len(dctMatArr) / inputDCTArrLen
	C.SVD(&dctMatArr[0], C.int(row), C.int(inputDCTArrLen))

}
