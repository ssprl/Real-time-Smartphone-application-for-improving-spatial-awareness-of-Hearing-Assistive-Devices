//
// Created by Abdullah Kucuk on 9/11/2017.
//


#include <cassert>
#include <cstring>
#include <jni.h>

#include <sys/types.h>
#include <SLES/OpenSLES.h>


#include "audio_main.h"
double ake_sum=0;
double ake_timing=0;

extern "C" {

JNIEXPORT jlong JNICALL
        Java_com_example_axk166230_TwoChanDOA_1v4_SignalProcessing_paramInitilization(JNIEnv *env, jobject thiz,int frequency, int stepSize,int frameSize, float stepFactor, int noisetype, float ThreadTime, float DurationTime, bool isEnchanced );
JNIEXPORT void JNICALL
        Java_com_example_axk166230_TwoChanDOA_1v4_SignalProcessing_realtimeProcessing(JNIEnv *env, jobject thiz,  jlong memoryPointer, jshortArray input);
JNIEXPORT void JNICALL
        Java_com_example_axk166230_TwoChanDOA_1v4_SignalProcessing_paramElimination(JNIEnv* env, jobject thiz, jlong memoryPointer);
JNIEXPORT jshortArray JNICALL
        Java_com_example_axk166230_TwoChanDOA_1v4_SignalProcessing_soundOutput(JNIEnv* env, jobject thiz, jlong memoryPointer, jint outputSelection);

JNIEXPORT jfloatArray JNICALL
        Java_com_example_axk166230_TwoChanDOA_1v4_SignalProcessing_angleOutput(JNIEnv* env, jobject thiz, jlong memoryPointer, jint outputSelection);

JNIEXPORT jshortArray JNICALL
        Java_com_example_axk166230_TwoChanDOA_1v4_SignalProcessing_dataOutput(JNIEnv* env, jobject thiz, jlong memoryPointer, jint outputSelection);
JNIEXPORT jfloat JNICALL
        Java_com_example_axk166230_TwoChanDOA_1v4_SignalProcessing_getTime(JNIEnv* env, jobject thiz, jlong memoryPointer);
JNIEXPORT jfloat JNICALL
        Java_com_example_axk166230_TwoChanDOA_1v4_SignalProcessing_getComputeTime(JNIEnv* env, jobject thiz, jlong memoryPointer);
JNIEXPORT jfloat JNICALL
        Java_com_example_axk166230_TwoChanDOA_1v4_SignalProcessing_getFilteringTime(JNIEnv* env, jobject thiz, jlong memoryPointer);






JNIEXPORT jboolean JNICALL
        Java_com_example_axk166230_twochandoa_1v4_MainActivity_createSLBufferQueueAudioPlayer(JNIEnv *env, jclass);
JNIEXPORT void JNICALL
        Java_com_example_axk166230_twochandoa_1v4_MainActivity_deleteSLBufferQueueAudioPlayer(JNIEnv *env, jclass type);

JNIEXPORT void JNICALL
        Java_com_example_axk166230_twochandoa_1v4_SignalProcessing_realtimeProcessing(JNIEnv *env, jobject thiz,  jlong memoryPointer, jshortArray input);
JNIEXPORT void JNICALL
        Java_com_example_axk166230_twochandoa_1v4_MainActivity_deleteAudioRecorder(JNIEnv *env, jclass type);
JNIEXPORT void JNICALL
        Java_com_example_axk166230_twochandoa_1v4_MainActivity_startPlay(JNIEnv *env, jclass type);
JNIEXPORT void JNICALL
        Java_com_example_axk166230_twochandoa_1v4_MainActivity_stopPlay(JNIEnv *env, jclass type);
}






JNIEXPORT jlong JNICALL
Java_com_example_axk166230_TwoChanDOA_1v4_SignalProcessing_paramInitilization(
        JNIEnv *env, jobject thiz,int frequency, int stepSize,int frameSize, float stepFactor, int noisetype, float ThreadTime, float DurationTime, bool isEnchanced){
    Variables *mainParam = (Variables*)malloc(sizeof(Variables));

    mainParam->stepSize = stepSize;
    mainParam->samplingRate = frequency;
    mainParam->frameSize = frameSize;  ////Threadshold Time
    mainParam->filLen = 1;
    mainParam->topMicBuffer = (float*)calloc(frameSize+mainParam->filLen-1,sizeof(float));
    mainParam->botMicBuffer = (float*)calloc(frameSize+mainParam->filLen-1,sizeof(float));

    mainParam->x_frame = (double*)malloc(stepSize * sizeof(double));
    mainParam->h_frame = (double*)malloc(stepSize * sizeof(double));

    mainParam->originalBuffer = (int*)calloc(2*stepSize,sizeof(int));
    mainParam->mixedBuffer = (short*)calloc(stepSize,sizeof(short));

    mainParam->Threadtime = ThreadTime;
    mainParam->DurationTime = DurationTime;
    mainParam->angle = 0;
    mainParam->angle_count = 0;
    mainParam->twoDOA = newDOA(DOA_NFFTT, stepSize, DOA_fss,ThreadTime,DurationTime);
    mainParam->isEnchanced=isEnchanced;
    mainParam->logMMSE=newparameters(DOA_NFFTT, stepSize, 50); //PERC=50
    mainParam->ake_counter=0;
    mainParam->ake_avg_timer=0;

    /*Memory Update*/
    mainParam->_output=(float*)calloc(3,sizeof(float)); // We have 3 output

    __android_log_print(ANDROID_LOG_ERROR, "Parameter Initialized","Successfully--IT WAS HEREEEE");
    return (jlong)mainParam;
}


JNIEXPORT void JNICALL
Java_com_example_axk166230_TwoChanDOA_1v4_SignalProcessing_realtimeProcessing(JNIEnv *env, jobject thiz,  jlong memoryPointer, jshortArray input){
    Variables *mainParam = (Variables*) memoryPointer;
    //startTimer(mainParam->timer);
    mainParam->_in = env->GetShortArrayElements(input, NULL);
    int i,j;//stepSize,frameSize,filLen;
    

    for (int k = 0; k < DOA_nn; k++)
    {

        mainParam-> x_frame[k] = mainParam->botMicBuffer[k];
        mainParam->h_frame[k] = mainParam->topMicBuffer[k];

    }
   //*
    clock_t t;
    t = clock();

    mainParam->twoDOA->doIT(mainParam->twoDOA,mainParam->x_frame, mainParam->h_frame, mainParam->angle_count, mainParam->twoDOA->prevmag_fft_Framex, mainParam->twoDOA->prevmag_fft_Frameh, mainParam->twoDOA->prevcorrtheta_est, mainParam->twoDOA->SFxmax, mainParam->twoDOA->SFxavg, mainParam->twoDOA->flagSFx, DOA_fss, DOA_nn, DOA_NFFTT,mainParam->Threadtime,mainParam->DurationTime,mainParam->logMMSE,mainParam->isEnchanced);

    t = clock() - t;

    ake_timing=(((double)t)/ CLOCKS_PER_SEC)+(mainParam->ake_avg_timer*(mainParam->ake_counter));

    mainParam->ake_counter= mainParam->ake_counter+1;
    mainParam->ake_avg_timer=ake_timing/(mainParam->ake_counter);//*/
     /* if (mainParam->angle_count %20 == 0){
    __android_log_print(ANDROID_LOG_INFO,"Time", "Time1= %1.9g",((float)mainParam->ake_avg_timer));
          }*/
    //__android_log_print(ANDROID_LOG_INFO,"Time", "Time1= %1.9g",(((double)t)/ CLOCKS_PER_SEC));


    /*if(mainParam->angle_count==100){

        tc=clock();
    }
*/
    /* FILE *file1 = fopen("/sdcard/mic.txt", "w+");

     fprintf(file1, "%d\t%ld click\t %1.9g\n", mainParam->angle_count, t,((float)t)/ CLOCKS_PER_SEC);
     fflush(file1);
     fclose(file1);*/



    //__android_log_print(ANDROID_LOG_INFO,"Angle", "The angle = %lf",mainParam->twoDOA->corrtheta_est);
    //__android_log_print(ANDROID_LOG_INFO,"Time", "Time= %1.9g",((float)t)/ CLOCKS_PER_SEC);
    mainParam->angle_count++;

    for(i = 0;i<mainParam->filLen-1;i++){
        mainParam->topMicBuffer[i] = mainParam->topMicBuffer[i+mainParam->stepSize];
        mainParam->botMicBuffer[i] = mainParam->botMicBuffer[i+mainParam->stepSize];
    }
    for(i = mainParam->filLen-1,j=0;i<mainParam->filLen+mainParam->stepSize-1;i++,j+=2){
        mainParam->topMicBuffer[i]= mainParam->topMicBuffer[i+mainParam->stepSize];
        mainParam->topMicBuffer[i+mainParam->stepSize]= mainParam->_in[j+1]*FAC;

        mainParam->botMicBuffer[i]= mainParam->botMicBuffer[i+mainParam->stepSize];
        mainParam->botMicBuffer[i+mainParam->stepSize]= mainParam->_in[j]*FAC;

        mainParam->originalBuffer[j] = mainParam->_in[j];
        mainParam->originalBuffer[j+1] = mainParam->_in[j+1];
    }
    for(i=0;i<mainParam->stepSize;i++){
        mainParam->mixedBuffer[i] = (mainParam->_in[i*2]+mainParam->_in[i*2+1])/2;
    }
    //__android_log_print(ANDROID_LOG_ERROR, "Parameter Computed 1st","Successfully");
    env->ReleaseShortArrayElements(input, mainParam->_in, 0);
    //tempEng = 0;

   
}


JNIEXPORT void JNICALL
Java_com_example_axk166230_TwoChanDOA_1v4_SignalProcessing_paramElimination(JNIEnv* env, jobject thiz, jlong memoryPointer){
    Variables *mainParam = (Variables*) memoryPointer;
    if(mainParam != NULL){
        /*
        if(mainParam->topMicBuffer!= NULL)
        free(mainParam->topMicBuffer);mainParam->topMicBuffer = NULL;
        if(mainParam->botMicBuffer!= NULL)
        free(mainParam->botMicBuffer);mainParam->botMicBuffer = NULL;
        //free(mainParam->outputBuffer);mainParam->outputBuffer = NULL;
        if(mainParam->x_frame!= NULL)
        free(mainParam->x_frame);mainParam->x_frame = NULL;
        if(mainParam->h_frame!= NULL)
        free(mainParam->h_frame);mainParam->h_frame = NULL;
        if(mainParam->originalBuffer!= NULL)
        free(mainParam->originalBuffer);mainParam->originalBuffer = NULL;
        if(mainParam->mixedBuffer!= NULL)
        free(mainParam->mixedBuffer);mainParam->mixedBuffer = NULL;
        if(mainParam->logMMSE!= NULL)
        free(mainParam->logMMSE);mainParam->logMMSE = NULL;
        if(mainParam->_output!= NULL)
        free(mainParam->_output);mainParam->_output = NULL;
        if(mainParam->_in!= NULL)
        free(mainParam->_in);mainParam->_in = NULL;
        //if(mainParam->twoDOA!= NULL)
        //free(mainParam->twoDOA);mainParam->twoDOA = NULL;
        //*/

        /*
        free(mainParam->y_prev);mainParam->y_prev = NULL;
        free(mainParam->y_curr);mainParam->y_curr = NULL;
        free(mainParam->y);mainParam->y = NULL;
        free(mainParam->e);mainParam->e = NULL;
        free(mainParam->w);mainParam->w = NULL;//*/

        //destroylogMMSE(&(mainParam->logmmsePtr));
        //destroyTimer(&(mainParam->timer));
        //free(mainParam);mainParam = NULL;
    }
}

JNIEXPORT jshortArray JNICALL
Java_com_example_axk166230_TwoChanDOA_1v4_SignalProcessing_soundOutput(JNIEnv* env, jobject thiz, jlong memoryPointer, jint outputSelection){
    Variables* mainParam = (Variables*) memoryPointer;
    jshortArray output = env->NewShortArray(mainParam->stepSize);
    short *_output = env->GetShortArrayElements( output, NULL);
    int i;
    switch (outputSelection){
        case 0:		// Original

            for (i=0;i<mainParam->stepSize;i++){
                _output[i] = (short)checkRange(mainParam->mixedBuffer[i]);
            }


            for (i=0;i<mainParam->stepSize;i++){
                _output[i] = mainParam->mixedBuffer[i];
            }
            break;

    }
    env->ReleaseShortArrayElements(output, _output, 0);
    return output;
}


JNIEXPORT jfloatArray JNICALL
Java_com_example_axk166230_TwoChanDOA_1v4_SignalProcessing_angleOutput(JNIEnv* env, jobject thiz, jlong memoryPointer, jint outputSelection){
    Variables* mainParam = (Variables*) memoryPointer;
    jfloatArray output = env->NewFloatArray(3);
    //float tmp = 0;
    mainParam->_output = env->GetFloatArrayElements( output, NULL);
    //int i;
    switch (outputSelection){
        case 0:		// Original



            //tmp = max1(2,3);

            //mainParam->angle = tmp;
            //mainParam->angle_count++;
            mainParam->_output[0] = mainParam->twoDOA->corrtheta_est;
            mainParam->_output[1] = mainParam->twoDOA->SFxavg;
            mainParam->_output[2] = mainParam->ake_counter;

            break;
        case 1:		// Enhanced
            mainParam->_output[0] = -1;
            /*
            for (i=0,j=0;i<2*mainParam->stepSize;i+=2,j++){
                _output[i] = (short)checkRange(mainParam->e[j]);
                _output[i+1] = (short)checkRange(mainParam->e[j]);
            }*/
            break;
        case 2:
            break;
    }
    env->ReleaseFloatArrayElements(output, mainParam->_output, 0);
//    free(_output);
  //  _output=NULL;
    return output;

}

JNIEXPORT jshortArray JNICALL
Java_com_example_axk166230_TwoChanDOA_1v4_SignalProcessing_dataOutput(JNIEnv* env, jobject thiz, jlong memoryPointer, jint outputSelection){
    Variables* mainParam = (Variables*) memoryPointer;
    jshortArray output = env->NewShortArray(2*mainParam->stepSize);
    short *_output = env->GetShortArrayElements(output, NULL);
    int i;
    switch (outputSelection){
        case 0:
            for(i= 0;i<2*mainParam->stepSize;i++){
                _output[i] = mainParam->originalBuffer[i];
            }
            break;
      /*  case 1:
            for(i= 0;i<2*mainParam->stepSize;i++){
                _output[i] = (short)checkRange(mainParam->outputBuffer[i]);
            }
            break;
        case 2:
            break;//*/
    }
    env->ReleaseShortArrayElements(output, _output, 0);
    return output;
}














int checkRange(float input)
{
    int output;
    if(input>1){
        output = 32767;
    }else if(input<-1){
        output = -32768;
    }else
        output = 32768*input;
    return output;
}
//////////////
///////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*JNIEXPORT jboolean JNICALL
Java_com_google_sample_echo_MainActivity_createSLBufferQueueAudioPlayer(JNIEnv *env, jclass type) {
    SampleFormat sampleFormat;
    memset(&sampleFormat, 0, sizeof(sampleFormat));
    sampleFormat.pcmFormat_ = (uint16_t)engine.bitsPerSample_;
    sampleFormat.framesPerBuf_ = engine.fastPathFramesPerBuf_;

    // SampleFormat.representation_ = SL_ANDROID_PCM_REPRESENTATION_SIGNED_INT;
    sampleFormat.channels_ = (uint16_t)engine.sampleChannels_;
    sampleFormat.sampleRate_ = engine.fastPathSampleRate_;

    engine.player_ = new AudioPlayer(&sampleFormat, engine.slEngineItf_);
    assert(engine.player_);
    if(engine.player_ == nullptr)
        return JNI_FALSE;

    engine.player_->SetBufQueue(engine.recBufQueue_, engine.freeBufQueue_);
    engine.player_->RegisterCallback(EngineService, (void*)&engine);

    return JNI_TRUE;
}

JNIEXPORT void JNICALL
Java_com_google_sample_echo_MainActivity_deleteSLBufferQueueAudioPlayer(JNIEnv *env, jclass type) {
    if(engine.player_) {
        delete engine.player_;
        engine.player_= nullptr;
    }
}



JNIEXPORT void JNICALL
Java_com_google_sample_echo_MainActivity_deleteAudioRecorder(JNIEnv *env, jclass type) {
    if(engine.recorder_)
        delete engine.recorder_;

    engine.recorder_ = nullptr;
}

JNIEXPORT void JNICALL
Java_com_google_sample_echo_MainActivity_startPlay(JNIEnv *env, jclass type) {

    engine.frameCount_  = 0;*/
    /*
     * start player: make it into waitForData state
     */
 /*
    if(SL_BOOLEAN_FALSE == engine.player_->Start()){
        LOGE("====%s failed", __FUNCTION__);
        return;
    }
    engine.recorder_->Start();
}

JNIEXPORT void JNICALL
Java_com_google_sample_echo_MainActivity_stopPlay(JNIEnv *env, jclass type) {
    engine.recorder_->Stop();
    engine.player_ ->Stop();

    delete engine.recorder_;
    delete engine.player_;
    engine.recorder_ = NULL;
    engine.player_ = NULL;
}



uint32_t dbgEngineGetBufCount(void) {
    uint32_t count = engine.player_->dbgGetDevBufCount();
    count += engine.recorder_->dbgGetDevBufCount();
    count += engine.freeBufQueue_->size();
    count += engine.recBufQueue_->size();

    LOGE("Buf Disrtibutions: PlayerDev=%d, RecDev=%d, FreeQ=%d, "
                 "RecQ=%d",
         engine.player_->dbgGetDevBufCount(),
         engine.recorder_->dbgGetDevBufCount(),
         engine.freeBufQueue_->size(),
         engine.recBufQueue_->size());
    if(count != engine.bufCount_) {
        LOGE("====Lost Bufs among the queue(supposed = %d, found = %d)",
             BUF_COUNT, count);
    }
    return count;
}*/

/*
 * simple message passing for player/recorder to communicate with engine
 */

 /*
bool EngineService(void* ctx, uint32_t msg, void* data ) {
    assert(ctx == &engine);
    switch (msg) {
        case ENGINE_SERVICE_MSG_KICKSTART_PLAYER:
            engine.player_->PlayAudioBuffers(PLAY_KICKSTART_BUFFER_COUNT);
            // we only allow it to call once, so tell caller do not call
            // anymore
            return false;
        case ENGINE_SERVICE_MSG_RETRIEVE_DUMP_BUFS:
            *(static_cast<uint32_t*>(data)) = dbgEngineGetBufCount();
            break;
        default:
            assert(false);
            return false;
    }

    return true;
}
*/


