package com.example.axk166230.TwoChanDOA_v4;

import android.app.Activity;
import android.util.Log;

import java.util.concurrent.BlockingQueue;

/**
 * Created by Abdullah Kucuk on 9/11/2017.
 */


public class Processing extends Activity implements Runnable {




    MainActivity myActivity = new MainActivity();


    private BlockingQueue<AudioFrame> input;

    private BlockingQueue<AudioFrame> output;
    private SignalProcessing signalProcessor;
    private Thread speechThread;
    private int counter;
    private float[] temp_angle;
    private float angle;
    private float SFxavg;
    private int ake_counter;

    private static String TAG = Processing.class.getSimpleName();
    public Processing(BlockingQueue<AudioFrame> input, BlockingQueue<AudioFrame> output){
        this.input = input;

        this.output = output;
        signalProcessor = new SignalProcessing();
        speechThread = new Thread(this);
        speechThread.start();
        Log.e(TAG, "onCreate initialize successfully!");
    }
    @Override
    public void run() {
        try{
            loop:while(true){
                AudioFrame currentFrame = null;



                currentFrame = input.take();



                if(currentFrame==Settings.STOP){
                  /*  Log.e(TAG,"Stage 1 compute time:"+signalProcessor.getFilteringTime()+" ms");
                    Log.e(TAG,"Stage 2 compute time:"+signalProcessor.getComputeTime()+" ms");*/
                   // signalProcessor.release();
                    output.put(currentFrame);
                    break loop;
                }
                //Log.e(TAG, "frame taken successfully!");

                signalProcessor.frameProcess(currentFrame.getAudio());
                //Log.e(TAG, "frame Processing successfully!");

                if(Settings.playback){
                    //currentFrame.setAduio(signalProcessor.getSoundOutput());
                    temp_angle = signalProcessor.getAngleOutput();
                    angle = temp_angle[0];
                    Settings.output_angle = angle;
                    SFxavg = temp_angle[1];
                    Settings.SFx_avg = SFxavg;
                    //angle = temp_angle[0];
                    ake_counter=(int)temp_angle[2];



                    if (ake_counter%5==0){
                        myActivity.runOnUiThread(new Runnable() {
                            public void run() {
                                myActivity.updateGUI();
                            }
                        });}//*/

                }
                /*
                if(Settings.debugLevel==0){
                    currentFrame.setDebug(signalProcessor.getDataOutput());
                }
                */
                //output.put(currentFrame);
            }
        }catch(InterruptedException e){
            Thread.currentThread().interrupt();
            e.printStackTrace();
        }
    }
}
