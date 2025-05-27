#include <stdio.h>
#include <iostream>
#include <math.h>
#include <cmath>

#define EMSCRIPTEN_KEEPALIVE __attribute__((used))
#define PI 3.14159265358979323846264338327950288

#define WEBAUDIO_BLOCK_SIZE 128
#define BUFFERED_BLOCK_SIZE 2048

unsigned int uleftshift(int v, int shift) {
   return static_cast<unsigned int>(v) << shift;
}

unsigned int urightshift(int v, int shift) {
   return static_cast<unsigned int>(v) >> shift;
}

class FFT {
   private:
   int size;
   int _csize;
   std::vector<float> table;
   std::vector<unsigned int> _bitrev;
   int _width;
   std::vector<float> _out;
   std::vector<float> _data;
   int _inv;

   public:
   // Constructor
   FFT(int size){
      this->size = size;
      if (this->size <= 1 || (this->size & (this->size - 1)) != 0)
         throw std::invalid_argument("FFT size must be a power of two and bigger than 1");
      
      this->_csize = size << 1;

      this->table.resize(this->size*2);
      for (int i = 0; i < (int)table.size(); i += 2){
         const float angle = PI * i / this->size;
         table[i] = cosf(angle);
         table[i + 1] = -sinf(angle);
      }
      //this->table = table;

      int power = 0;
      for (int t = 1; this->size > t; t <<= 1)
         power++;

      this->_width = power % 2 == 0 ? power - 1 : power;
      
      this->_bitrev.resize(1 << this->_width);
      for (int j = 0; j < (int)this->_bitrev.size(); j++){
         this->_bitrev[j] = 0;
         for (int shift = 0; shift < this->_width; shift += 2){
            int revShift = this->_width - shift - 2;
            this->_bitrev[j] |= ( urightshift(j, shift) & 3 ) << revShift;
         }
      }

      this->_out = std::vector<float>();
      this->_data = std::vector<float>();
      this->_inv = 0;
   }

   // Destructor
   ~FFT() {
      std::cout << "FFT finished" << std::endl;
   }

   void fromComplexArray(const std::vector<float>& complex, std::vector<float>& storage) {
      if (storage.empty()) {
         storage.resize(complex.size() / 2);
      }

      for (int i = 0; i < (int)complex.size(); i += 2) {
         storage[i / 2] = complex[i];
      }
   }

   std::vector<float> createComplexArray() {
      return std::vector<float>(this->_csize, 0.0f);
   }

   std::vector<float> toComplexArray(const std::vector<float>& input, std::vector<float>& storage){
      std::vector<float> result;

      if (!storage.empty()){
         result = storage;
      }else result = this->createComplexArray();

      for (int i = 0; i < (int)result.size(); i += 2){
         result[i] = input[urightshift(i, 1)];
         result[i+1] = 0;
      }

      return result;
   }

   void completeSpectrum(std::vector<float>& spectrum){
      int size = this->_csize;
      unsigned int half = urightshift(size, 1);
      for (unsigned int i = 2; i < half; i += 2){
         spectrum[size - i] = spectrum[i];
         spectrum[size - i + 1] = -spectrum[i + 1];
      }
   }

   void transform(std::vector<float>& out, std::vector<float>& data){
      //if (out.size() == data.size() && std::equal(out.begin(), out.end(), data.begin()))
      //   throw std::invalid_argument("Input and output buffers must be different");

      this->_out = out;
      this->_data = data;
      this->_inv = 0;
      this->_transform4();
      out = this->_out;
      this->_out.clear();
      this->_data.clear();
   }

   void realTransform(std::vector<float>& out, std::vector<float>& data){
      //if (out.size() == data.size() && std::equal(out.begin(), out.end(), data.begin()))
      //   throw std::invalid_argument("Input and output buffers must be different");

      // tengo que revisar porque fft no se aplica bien

      this->_out = out;
      this->_data = data;
      this->_inv = 0;
      this->_realTransform4();
      out = this->_out; // es necesario
      //this->_out.clear();
      //this->_data.clear();
   }

   void inverseTransform(std::vector<float>& out, std::vector<float>& data){
      //if (out.size() == data.size() && std::equal(out.begin(), out.end(), data.begin()))
      //   throw std::invalid_argument("Input and output buffers must be different");

      this->_out = out;
      this->_data = data;
      this->_inv = 1;
      this->_transform4();
      out = this->_out; // es necesario
      for (int i = 0; i < (int)out.size(); i++)
         out[i] /= this->size;
      //this->_out.clear();
      //this->_data.clear();
   }

   void _transform4(){
      std::vector<float>& out = this->_out;
      int size = this->_csize;

      int width = this->_width;
      int step = 1 << width;
      int len = (size / step) << 1;

      int outOff;
      int t;
      std::vector<unsigned int>& bitrev = this->_bitrev;
      if (len == 4){
         for (outOff = 0, t = 0; outOff < size; outOff += len, t++){
            this->_singleTransform2(outOff, bitrev[t], step);
         }
      }else{
         // len == 8
         for (outOff = 0, t = 0; outOff < size; outOff += len, t++){
            this->_singleTransform4(outOff, bitrev[t], step);
         }
      }

      const int inv = this->_inv ? -1 : 1;
      std::vector<float>& table = this->table;
      for (step >>= 2; step >= 2; step >>= 2){
         len = (size / step) << 1;
         unsigned int quarterLen = urightshift(len, 2);
         
         for (outOff = 0; outOff < size; outOff += len){
            int limit = outOff + quarterLen;
            for (int i = outOff, k = 0; i < limit; i += 2, k += step){
               const int A = i;
               const int B = A + quarterLen;
               const int C = B + quarterLen;
               const int D = C + quarterLen;

               // Original values
               const float Ar = out[A];
               const float Ai = out[A+1];
               const float Br = out[B];
               const float Bi = out[B+1];
               const float Cr = out[C];
               const float Ci = out[C+1];
               const float Dr = out[D];
               const float Di = out[D+1];

               // Middle values
               const float MAr = Ar;
               const float MAi = Ai;

               const float tableBr = table[k];
               const float tableBi = inv*table[k+1];
               const float MBr = Br*tableBr-Bi*tableBi;
               const float MBi = Br*tableBi+Bi*tableBr;

               const float tableCr = table[2*k];
               const float tableCi = inv*table[2*k+1];
               const float MCr = Cr*tableCr-Ci*tableCi;
               const float MCi = Cr*tableCi+Ci*tableCr;

               const float tableDr = table[3*k];
               const float tableDi = inv*table[3*k+1];
               const float MDr = Dr*tableDr-Di*tableDi;
               const float MDi = Dr*tableDi+Di*tableDr;

               // Pre-Final values
               const float T0r = MAr+MCr;
               const float T0i = MAi+MCi;
               const float T1r = MAr-MCr;
               const float T1i = MAi-MCi;
               const float T2r = MBr+MDr;
               const float T2i = MBi+MDi;
               const float T3r = inv*(MBr-MDr);
               const float T3i = inv*(MBi-MDi);

               // Final values
               const float FAr = T0r+T2r;
               const float FAi = T0i+T2i;

               const float FCr = T0r-T2r;
               const float FCi = T0i-T2i;

               const float FBr = T1r+T3i;
               const float FBi = T1i-T3r;

               const float FDr = T1r-T3i;
               const float FDi = T1i+T3r;

               out[A] = FAr;
               out[A+1] = FAi;
               out[B] = FBr;
               out[B+1] = FBi;
               out[C] = FCr;
               out[C+1] = FCi;
               out[D] = FDr;
               out[D+1] = FDi;
            }
         }
      }


   }

   void _realTransform4(){
      std::vector<float>& out = this->_out;
      int size = this->_csize;

      int width = this->_width;
      int step = 1 << width;
      int len = (size / step) << 1;

      int outOff;
      int t;
      std::vector<unsigned int>& bitrev = this->_bitrev;
      if (len == 4){
         for (outOff = 0, t = 0; outOff < size; outOff += len, t++){
            this->_singleRealTransform2(outOff, urightshift(bitrev[t], 1), urightshift(step, 1));
         }
      }else{
         // len == 8
         for (outOff = 0, t = 0; outOff < size; outOff += len, t++){
            this->_singleRealTransform4(outOff, urightshift(bitrev[t], 1), urightshift(step, 1));
         }
      }

      const int inv = this->_inv ? -1 : 1;
      std::vector<float>& table = this->table;
      for (step >>= 2; step >= 2; step >>= 2){
         len = (size / step) << 1;
         unsigned int halfLen = urightshift(len, 1);
         unsigned int quarterLen = urightshift(halfLen, 1);
         unsigned hquarterLen = urightshift(quarterLen, 1);

         for (outOff = 0; outOff < size; outOff += len){
            for (unsigned int i = 0, k = 0; i <= hquarterLen; i += 2, k += step){
               const int A = outOff + i;
               const int B = A + quarterLen;
               const int C = B + quarterLen;
               const int D = C + quarterLen;

               // Original values
               const float Ar = out[A];
               const float Ai = out[A+1];
               const float Br = out[B];
               const float Bi = out[B+1];
               const float Cr = out[C];
               const float Ci = out[C+1];
               const float Dr = out[D];
               const float Di = out[D+1];

               // Middle values
               const float MAr = Ar;
               const float MAi = Ai;

               const float tableBr = table[k];
               const float tableBi = inv*table[k+1];
               const float MBr = Br*tableBr-Bi*tableBi;
               const float MBi = Br*tableBi+Bi*tableBr;

               const float tableCr = table[2*k];
               const float tableCi = inv*table[2*k+1];
               const float MCr = Cr*tableCr-Ci*tableCi;
               const float MCi = Cr*tableCi+Ci*tableCr;

               const float tableDr = table[3*k];
               const float tableDi = inv*table[3*k+1];
               const float MDr = Dr*tableDr-Di*tableDi;
               const float MDi = Dr*tableDi+Di*tableDr;

               // Pre-Final values
               const float T0r = MAr+MCr;
               const float T0i = MAi+MCi;
               const float T1r = MAr-MCr;
               const float T1i = MAi-MCi;
               const float T2r = MBr+MDr;
               const float T2i = MBi+MDi;
               const float T3r = inv*(MBr-MDr);
               const float T3i = inv*(MBi-MDi);

               // Final values
               const float FAr = T0r+T2r;
               const float FAi = T0i+T2i;

               const float FBr = T1r+T3i;
               const float FBi = T1i-T3r;

               out[A] = FAr;
               out[A+1] = FAi;
               out[B] = FBr;
               out[B+1] = FBi;

               // Output final middle point
               if (i==0){
                  const float FCr = T0r-T2r;
                  const float FCi = T0i-T2i;
                  out[C] = FCr;
                  out[C + 1] = FCi;
                  continue;
               }

               // Do not overwrite ourselves
               if (i==hquarterLen)
                  continue;

               // In the flipped case:
               // MAi = -MAi
               // MBr=-MBi, MBi=-MBr
               // MCr=-MCr
               // MDr=MDi, MDi=MDr
               const float ST0r = T1r;
               const float ST0i = -T1i;
               const float ST1r = T0r;
               const float ST1i = -T0i;
               const float ST2r = -inv*T3i;
               const float ST2i = -inv*T3r;
               const float ST3r = -inv*T2i;
               const float ST3i = -inv*T2r;

               const float SFAr = ST0r+ST2r;
               const float SFAi = ST0i+ST2i;

               const float SFBr = ST1r+ST3i;
               const float SFBi = ST1i-ST3r;

               const float SA = outOff+quarterLen-i;
               const float SB = outOff+halfLen-i;

               out[SA] = SFAr;
               out[SA+1] = SFAi;
               out[SB] = SFBr;
               out[SB+1] = SFBi;
            }
         }
      }
   }

   void _singleTransform2(int outOff, int off, int step){
      std::vector<float>& out = this->_out;
      std::vector<float>& data = this->_data;

      const float evenR = data[off];
      const float evenI = data[off+1];
      const float oddR = data[off+step];
      const float oddI = data[off+step+1];

      const float leftR = evenR + oddR;
      const float leftI = evenI + oddI;
      const float rightR = evenR - oddR;
      const float rightI = evenI - oddI;

      out[outOff] = leftR;
      out[outOff+1] = leftI;
      out[outOff+2] = rightR;
      out[outOff+3] = rightI;
   }

   void _singleRealTransform2(int outOff, int off, int step){
      std::vector<float>& out = this->_out;
      std::vector<float>& data = this->_data;

      const float evenR = data[off];
      const float oddR = data[off+step];

      const float leftR = evenR + oddR;
      const float rightR = evenR - oddR;

      out[outOff] = leftR;
      out[outOff+1] = 0;
      out[outOff+2] = rightR;
      out[outOff+3] = 0;
   }

   void _singleTransform4(int outOff, int off, int step){
      std::vector<float>& out = this->_out;
      std::vector<float>& data = this->_data;
      const int inv = this->_inv ? -1 : 1;
      const int step2 = step*2;
      const int step3 = step*3;

      // Original values
      const float Ar = data[off];
      const float Ai = data[off+1];
      const float Br = data[off+step];
      const float Bi = data[off+step+1];
      const float Cr = data[off+step2];
      const float Ci = data[off+step2+1];
      const float Dr = data[off+step3];
      const float Di = data[off+step3+1];

      // Pre-Final values
      const float T0r = Ar+Cr;
      const float T0i = Ai+Ci;
      const float T1r = Ar-Cr;
      const float T1i = Ai-Ci;
      const float T2r = Br+Dr;
      const float T2i = Bi+Di;
      const float T3r = inv*(Br-Dr);
      const float T3i = inv*(Bi-Di);

      // Final values
      const float FAr = T0r+T2r;
      const float FAi = T0i+T2i;

      const float FBr = T1r+T3i;
      const float FBi = T1i-T3r;

      const float FCr = T0r-T2r;
      const float FCi = T0i-T2i;

      const float FDr = T1r-T3i;
      const float FDi = T1i+T3r;

      out[outOff] = FAr;
      out[outOff+1] = FAi;
      out[outOff+2] = FBr;
      out[outOff+3] = FBi;
      out[outOff+4] = FCr;
      out[outOff+5] = FCi;
      out[outOff+6] = FDr;
      out[outOff+7] = FDi;
   }

   void _singleRealTransform4(int outOff, int off, int step){
      std::vector<float>& out = this->_out;
      std::vector<float>& data = this->_data;
      const int inv = this->_inv ? -1 : 1;
      const int step2 = step*2;
      const int step3 = step*3;

      // Original values
      const float Ar = data[off];
      const float Br = data[off+step];
      const float Cr = data[off+step2];
      const float Dr = data[off+step3];

      // Pre-Final values
      const float T0r = Ar+Cr;
      const float T1r = Ar-Cr;
      const float T2r = Br+Dr;
      const float T3r = inv*(Br-Dr);

      // Final values
      const float FAr = T0r+T2r;
      const float FBr = T1r;
      const float FBi = -T3r;
      const float FCr = T0r-T2r;
      const float FDr = T1r;
      const float FDi = T3r;

      out[outOff] = FAr;
      out[outOff+1] = 0;
      out[outOff+2] = FBr;
      out[outOff+3] = FBi;
      out[outOff+4] = FCr;
      out[outOff+5] = 0;
      out[outOff+6] = FDr;
      out[outOff+7] = FDi;
   }


};

class PhaseVocoder {
   private:
   int nbInputs;
   int nbOutputs;
   int blockSize;
   int hopSize;
   float nbOverlaps;
   std::vector<std::vector<std::vector<float > > > inputBuffers;
   std::vector<std::vector<float* > > inputBuffersHead;
   std::vector<std::vector<std::vector<float > > > inputBuffersToSend;
   std::vector<std::vector<std::vector<float > > > outputBuffers;
   std::vector<std::vector<std::vector<float > > > outputBuffersToRetrieve;
   int fftSize;
   float timeCursor;
   std::vector<float> hannWindow;
   FFT* fft;
   std::vector<float> freqComplexBuffer;
   std::vector<float> freqComplexBufferShifted;
   std::vector<float> timeComplexBuffer;
   std::vector<float> magnitudes;
   std::vector<int> peakIndexes;
   int nbPeaks;

   std::vector<float> genHannWindow(int length) {
      std::vector<float> win(length);
      for (int i = 0; i < length; i++) {
         win[i] = 0.5f * (1 - std::cos(2 * PI * i / static_cast<float>(length)));
      }
      return win;
   }

   void freeHannWindow(float* win) {
      delete[] win;
   }

   void reallocateChannelsIfNeeded(std::vector<std::vector<std::vector<float > > >& inputs, std::vector<std::vector<std::vector<float > > >& outputs){
      for (int i=0;i<this->nbInputs;i++){
         int nbChannels = (int)inputs[i].size();
         if (nbChannels != (int)this->inputBuffers[i].size()){
            this->allocateInputChannels(i, nbChannels);
         }
      }

      for (int i=0;i<this->nbInputs;i++){
         int nbChannels = (int)outputs[i].size();
         if (nbChannels != (int)this->outputBuffers[i].size()){
            this->allocateOutputChannels(i, nbChannels);
         }
      }
   }

   void allocateInputChannels(int inputIndex, int nbChannels){
      this->inputBuffers.resize(this->nbInputs);
      this->inputBuffers[inputIndex].resize(nbChannels);
      for (int i = 0; i < nbChannels; i++) {
         this->inputBuffers[inputIndex][i].resize(this->blockSize + hopSize, 0.0f);
      }

      this->inputBuffersHead.resize(this->nbInputs);
      this->inputBuffersToSend.resize(this->nbInputs);

      this->inputBuffersHead[inputIndex].resize(nbChannels);
      this->inputBuffersToSend[inputIndex].resize(nbChannels);
      
      for (int i = 0; i < nbChannels; i++) {
         this->inputBuffersHead[inputIndex][i] = this->inputBuffers[inputIndex][i].data();
         this->inputBuffersToSend[inputIndex][i].resize(this->blockSize, 0.0f);
      }
   }

   void allocateOutputChannels(int outputIndex, int nbChannels){
      this->outputBuffers[outputIndex].resize(nbChannels);
      for (int i = 0; i < nbChannels; i++) {
         this->outputBuffers[outputIndex][i].resize(this->blockSize, 0.0f);
      }

      this->outputBuffersToRetrieve[outputIndex].resize(nbChannels);
      for (int i = 0; i < nbChannels; i++) {
         this->outputBuffersToRetrieve[outputIndex][i].resize(this->blockSize, 0.0f);
      }
   }

   void readInputs(std::vector<std::vector<std::vector<float > > >& inputs){
      // when playback is paused, we may stop receiving new samples
      if (inputs[0].size()==0 && inputs[0][0].size()==0){
         for (int i = 0; i < this->nbInputs; i++){
            for (int j = 0; j < (int)this->inputBuffers[i].size(); j++){
               this->inputBuffers[i][j].assign(this->blockSize, 0);
            }
         }
         return;
      }

      for (int i = 0; i < this->nbInputs; i++) {
         for (int j = 0; j < (int)this->inputBuffers[i].size(); j++) {
            const std::vector<float>& webAudioBlock = inputs[i][j];
            std::copy(
                  webAudioBlock.begin(), 
                  webAudioBlock.end(), 
                  this->inputBuffers[i][j].begin() + this->blockSize
            );
         }
      }
      
   }

   void writeOutputs(std::vector<std::vector<std::vector<float > > >& outputs){
      for (int i = 0; i < this->nbInputs; i++) {
         for (int j = 0; j < (int)this->inputBuffers[i].size(); j++) {
            std::vector<float> webAudioBlock(
               this->outputBuffers[i][j].begin(),
               this->outputBuffers[i][j].begin() + hopSize
            );

            outputs[i][j] = webAudioBlock;
         }
      }
   }

   void shiftInputBuffers(){
      for (int i = 0; i < this->nbInputs; i++) {
         for (int j = 0; j < (int)this->inputBuffers[i].size(); j++) {
            
            std::vector<float>& buf = this->inputBuffers[i][j];
            std::copy(
               buf.begin() + this->hopSize,
               buf.end(),
               buf.begin()
            );

            /*std::fill(
               buf.end() - this->hopSize,
               buf.end(),
               0.0f
            );*/

         }
      }
   }

   void shiftOutputBuffers(){
      for (int i = 0; i < this->nbOutputs; i++) {
         for (int j = 0; j < (int)this->outputBuffers[i].size(); j++) {
            std::vector<float>& buf = this->outputBuffers[i][j];

            std::copy(
                  buf.begin() + hopSize,
                  buf.end(),
                  buf.begin()
            );

            int fill_start = this->blockSize - hopSize;
            std::fill(
                  buf.begin() + fill_start,
                  buf.end(),
                  0.0f
            );
         }
      }
   }

   void prepareInputBuffersToSend(){
      for (int i = 0; i < this->nbInputs; i++){
         for (int j = 0; j < (int)this->inputBuffers[i].size(); j++){
            std::copy(
               this->inputBuffersHead[i][j],
               this->inputBuffersHead[i][j] + this->blockSize,
               this->inputBuffersToSend[i][j].begin()
            );
         }
      }
   }

   void handleOutputBuffersToRetrieve(){
      for (int i = 0; i < this->nbOutputs; i++){
         for (int j = 0; j < (int)this->outputBuffers[i].size(); j++){
            for (int k = 0; k < this->blockSize; k++){
               this->outputBuffers[i][j][k] += this->outputBuffersToRetrieve[i][j][k] / this->nbOverlaps;
            }
         }
      }
   }

   void applyHannWindow(std::vector<float>& input){
      for (int i = 0; i < this->blockSize; i++){
         input[i] = input[i] * this->hannWindow[i];
      }
   }

   void computeMagnitudes(){
      int i = 0, j = 0;
      while (i < (int)this->magnitudes.size()){
         float real = this->freqComplexBuffer[i];
         float imag = this->freqComplexBuffer[j + 1];
         this->magnitudes[i] = real * real + imag * imag;
         i+=1;
         j+=2;
      }
   }

   void findPeaks(){
      this->nbPeaks = 0;
      int i = 2;
      int end = (int)this->magnitudes.size()-2;

      while (i < end){
         float mag = this->magnitudes[i];

         if (this->magnitudes[i-1] >= mag || this->magnitudes[i-2] >= mag){
            i++;
            continue;
         }

         if (this->magnitudes[i+1] >= mag || this->magnitudes[i + 2] >= mag){
            i++;
            continue;
         }

         this->peakIndexes[this->nbPeaks] = i;
         this->nbPeaks++;
         i += 2;
      }
   }

   void shiftPeaks(float pitchFactor){
      this->freqComplexBufferShifted.assign(this->freqComplexBufferShifted.size(), 0.0f);

      for (int i=0; i < this->nbPeaks; i++){
         int peakIndex = this->peakIndexes[i];
         int peakIndexShifted = roundf(peakIndex * pitchFactor);

         if (peakIndexShifted > (int)this->magnitudes.size()){
            break;
         }

         int startIndex = 0;
         int endIndex = this->fftSize;
         if (i > 0){
            int peakIndexBefore = this->peakIndexes[i-1];
            startIndex = peakIndex - floorf((peakIndex - peakIndexBefore) / 2.0f);
         }
         if (i < this->nbPeaks - 1){
            int peakIndexAfter = this->peakIndexes[i+1];
            endIndex = peakIndex + ceilf((peakIndexAfter - peakIndex) / 2.0f);
         }

         int startOffset = startIndex - peakIndex;
         int endOffset = endIndex - peakIndex;

         for (int j = startOffset; j < endOffset; j++){
            int binIndex = peakIndex + j;
            int binIndexShifted = peakIndexShifted + j;
            
            if (binIndexShifted >= (int)this->magnitudes.size()){
               break;
            }

            float omegaDelta = 2 * PI * (binIndexShifted - binIndex) / this->fftSize;
            float phaseShiftReal = cosf(omegaDelta * this->timeCursor);
            float phaseShiftImag = sinf(omegaDelta * this->timeCursor);

            int indexReal = binIndex * 2;
            int indexImag = indexReal + 1;
            float valueReal = this->freqComplexBuffer[indexReal];
            float valueImag = this->freqComplexBuffer[indexImag];

            float valueShiftedReal = valueReal * phaseShiftReal - valueImag * phaseShiftImag;
            float valueShiftedImag = valueReal * phaseShiftImag + valueImag * phaseShiftReal;

            int indexShiftedReal = binIndexShifted * 2;
            int indexShiftedImag = indexShiftedReal + 1;

            this->freqComplexBufferShifted[indexShiftedReal] += valueShiftedReal;
            this->freqComplexBufferShifted[indexShiftedImag] += valueShiftedImag;
         }

      }
   }

   public:
   // Constructor
   PhaseVocoder(int ninputs, int noutputs){
      this->nbInputs = ninputs;
      this->nbOutputs = noutputs;

      this->blockSize = BUFFERED_BLOCK_SIZE;

      this->hopSize = WEBAUDIO_BLOCK_SIZE;

      this->nbOverlaps = this->blockSize / this->hopSize;

      this->inputBuffers.resize(this->nbInputs);
      this->inputBuffersHead.resize(this->nbInputs);
      this->inputBuffersToSend.resize(this->nbInputs);
      for (int i=0;i<this->nbInputs;i++){
         this->allocateInputChannels(i, 1);
      }

      this->outputBuffers.resize(this->nbInputs);
      this->outputBuffersToRetrieve.resize(this->nbInputs);
      for (int i=0;i<this->nbOutputs;i++){
         this->allocateOutputChannels(i, 1);
      }

      this->fftSize = this->blockSize;

      this->timeCursor = 0;

      this->hannWindow = this->genHannWindow(this->blockSize);

      this->fft = new FFT(this->fftSize);
      this->freqComplexBuffer = this->fft->createComplexArray();
      this->freqComplexBufferShifted = this->fft->createComplexArray();
      this->timeComplexBuffer = this->fft->createComplexArray();
      this->magnitudes.resize(this->fftSize / 2 + 1, 0.0f);
      this->peakIndexes.resize(this->magnitudes.size());
      this->nbPeaks = 0;
      
   }

   // Destructor
   ~PhaseVocoder() {
      //delete[] hannWindow;
      delete fft;
      std::cout << "Phase Vocoder finished" << std::endl;
   }

   void processOLA(std::vector<std::vector<std::vector<float > > >& inputs, std::vector<std::vector<std::vector<float > > >& outputs, float pitch){
      
      for (int i=0;i<this->nbInputs;i++){
         for (int j=0;j<(int)inputs[i].size();j++){
            std::vector<float>& input = inputs[i][j];
            std::vector<float>& output = outputs[i][j];

            this->applyHannWindow(input);

            this->fft->realTransform(this->freqComplexBuffer, input);

            this->computeMagnitudes();
            this->findPeaks();
            this->shiftPeaks(pitch);

            this->fft->completeSpectrum(this->freqComplexBufferShifted);
            this->fft->inverseTransform(this->timeComplexBuffer, this->freqComplexBufferShifted);
            this->fft->fromComplexArray(this->timeComplexBuffer, output);

            this->applyHannWindow(output);
         }
      }

      this->timeCursor += this->hopSize;

   }

   void process(float* inputs, float* outputs, int nbInput, int channelCount, int frameCount, float pitch){
      std::vector<std::vector<std::vector<float > > > input;
      input.resize(nbInput);
      std::vector<std::vector<std::vector<float > > > output;
      output.resize(nbInput);

      const int samplesPerInput = channelCount * frameCount;

      for (int n=0;n<nbInput;n++){
         input[n].resize(channelCount);
         output[n].resize(channelCount);
         const int inputOffset = n * samplesPerInput;

         for (int channel = 0; channel < channelCount; channel++){
            input[n][channel].resize(frameCount);
            output[n][channel].resize(frameCount);
            const int channelOffset = inputOffset + channel * frameCount;

            for (int i=0;i<frameCount;i++){
               input[n][channel][i] = inputs[channelOffset + i];
               output[n][channel][i] = outputs[channelOffset + i];
            }
         }
      }

      this->reallocateChannelsIfNeeded(input, output);

      this->readInputs(input);
      this->shiftInputBuffers();
      this->prepareInputBuffersToSend();
      this->processOLA(this->inputBuffersToSend, this->outputBuffersToRetrieve, pitch);
      this->handleOutputBuffersToRetrieve();
      this->writeOutputs(output);
      this->shiftOutputBuffers();

      // Volver a pasar output[n][channel][frame] a outputs array 1D
      for (int n = 0; n < nbInput; n++){
         for (int c = 0; c < channelCount; c++){
            const int startIdx = n * samplesPerInput + c * frameCount;
            for (int i = 0; i < frameCount; i++){
               outputs[startIdx+i] = output[n][c][i];
            }
         }
      }

   }

};

PhaseVocoder* pv;

extern "C" {
   #ifndef EMSCRIPTEN_KEEPALIVE
   EMSCRIPTEN_KEEPALIVE
   #endif
   void create_controller(int niput, int noutput){
      pv = new PhaseVocoder(niput, noutput);
   }

   #ifndef EMSCRIPTEN_KEEPALIVE
   EMSCRIPTEN_KEEPALIVE
   #endif
   void destroy_controller(){
      delete pv;
   }

   #ifndef EMSCRIPTEN_KEEPALIVE
   EMSCRIPTEN_KEEPALIVE
   #endif
   void process(float* input, float* output, int nb, int channelCount, int frameCount, float pitch) {
      // input y output vienen en forma de un array 1D
      if (pv == nullptr || pitch == 1.0f){
         // No es necesario ajustar, copiar directamente
         memcpy(output, input, channelCount * frameCount * sizeof(float));
      }else{
         pv->process(input, output, nb, channelCount, frameCount, pitch);
      }

   }

}
