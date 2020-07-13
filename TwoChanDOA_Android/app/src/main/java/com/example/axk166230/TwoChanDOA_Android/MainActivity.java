package com.example.axk166230.TwoChanDOA_v4;

/*
 * Copyright 2015 The Android Open Source Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

import android.Manifest;
import android.content.Context;
import android.content.pm.PackageManager;
import android.graphics.Color;
import android.os.Bundle;
import android.support.v4.app.ActivityCompat;
import android.support.v7.app.AppCompatActivity;
import android.util.Log;
import android.view.View;
import android.widget.Button;
import android.widget.EditText;
import android.widget.ImageButton;
import android.widget.RadioButton;
import android.widget.RadioGroup;
import android.widget.TextView;

import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.atomic.AtomicBoolean;

import static android.content.ContentValues.TAG;





public class MainActivity extends AppCompatActivity
        implements ActivityCompat.OnRequestPermissionsResultCallback {
    private int mode;
    private String dialogFile;
    private AtomicBoolean processingFlag;
    private RadioGroup myradiogroup1;
    private RadioGroup myradiogroup2;
    private RadioButton rbEnhanced;
    private RadioButton rbOriginal;

    private RadioButton rbBabble;
    private RadioButton rbDrivingcar;
    private RadioButton rbMachinery;
    private EditText etFrameTime;
    public EditText etAngle;
    public EditText etSFxavg;
    public EditText etThreadDuration;
    //private RadioButton myradiobutton;
    private WaveRecorder waveRecorder;

   // private AudioReader audioReader;
    private BlockingQueue<AudioFrame> inputFrames;
    // private BlockingQueue<AudioFrame> inputFrames2;
    private BlockingQueue<AudioFrame> processedFrames;
    private int noisetype;
    private int playback;
    private static Context mContext;
    public Processing ake;
    private boolean isEnchanced=false;
    public TextView ake_orientation;
    private SignalProcessing signalProcessing;


    private static final int AUDIO_ECHO_REQUEST = 0;


    //private static MainActivity ins;


    public static Context getContext() {
        Log.i(TAG, "getContext: " + mContext);
        return mContext;
    }

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
       // getSupportActionBar().hide(); // for hiding title

        setContentView(R.layout.activity_main);


        Settings.setCallbackInterface(this);
        //Utilities.prepareDirectory(getString(R.string.appDirectory));

        inputFrames = new ArrayBlockingQueue<AudioFrame>(Settings.queueSize);
        // inputFrames2 = new ArrayBlockingQueue<AudioFrame>(Settings.queueSize);
        processedFrames = new ArrayBlockingQueue<AudioFrame>(Settings.queueSize);
        processingFlag = new AtomicBoolean();
        etFrameTime = (EditText)findViewById(R.id.et_frametime);
        etAngle = (EditText)findViewById(R.id.et_angle);
        etSFxavg = (EditText)findViewById(R.id.et_SFxavg);
        etThreadDuration = (EditText)findViewById(R.id.et_threadduration);
        Circle circle = (Circle) findViewById(R.id.circle);
        ake_orientation=(TextView)findViewById(R.id.textViewOrientation) ;
        signalProcessing= new SignalProcessing();




        enableButtons(false);
        enbaleSettings(false);
        setupListeners();

        Settings.et_degree = etAngle;
        Settings.et_SFxavg = etSFxavg;
        Settings.et_circle = circle;
        Settings.ake_orientation=ake_orientation;

        updateGUI();
    }

    private void enableButtons(boolean flag) {
        ((ImageButton)findViewById(R.id.ib_mic)).setEnabled(!flag);
        ((ImageButton)findViewById(R.id.ib_stop)).setEnabled(flag);
        ((ImageButton)findViewById(R.id.ib_settings)).setEnabled(!flag);
        ((Button)findViewById(R.id.b_enchanced)).setEnabled(!flag);
    }

    private void setupListeners(){
        ((ImageButton)findViewById(R.id.ib_mic)).setOnClickListener(buttonClick);
        ((ImageButton)findViewById(R.id.ib_stop)).setOnClickListener(buttonClick);
        ((ImageButton)findViewById(R.id.ib_settings)).setOnClickListener(buttonClick);
        ((Button)findViewById(R.id.b_save)).setOnClickListener(buttonClick);
        ((Button)findViewById(R.id.b_reset)).setOnClickListener(buttonClick);
        ((Button)findViewById(R.id.b_enchanced)).setOnClickListener(buttonClick);
    }

    private void enbaleSettings(boolean flag) {
        ((EditText) findViewById(R.id.et_frametime)).setEnabled(flag);
        ((EditText) findViewById(R.id.et_angle)).setEnabled(flag);
        ((EditText) findViewById(R.id.et_SFxavg)).setEnabled(flag);
        ((EditText) findViewById(R.id.et_threadduration)).setEnabled(flag);
        ((Button) findViewById(R.id.b_save)).setEnabled(flag);
        ((Button) findViewById(R.id.b_reset)).setEnabled(flag);

    }

    private View.OnClickListener buttonClick = new View.OnClickListener(){
        public void onClick(View v){
            switch (v.getId()){
                case R.id.ib_mic:{
                    //tv.setText(stringFromJNI());
                    StartRecord();


                    ////////////////////////
                    //////////////////////
                    Circle circle = (Circle) findViewById(R.id.circle);
                    //CircleView circle=new CircleView(this,45) ;
                    circleanimation animation = new circleanimation(circle, Settings.output_angle);
                    //animation.setInterpolator(new BounceInterpolator());
                    animation.setDuration(1000);
                    circle.startAnimation(animation);

                    String output_a = Float.toString(Settings.output_angle);
                    etAngle.setText(output_a);
                    etSFxavg.setText("0");                  ////////////////////////////
                    ///////////////////////////


                    break;
                }
                case R.id.ib_settings:{
                    ((ImageButton)findViewById(R.id.ib_mic)).setEnabled(false);
                    ((ImageButton)findViewById(R.id.ib_stop)).setEnabled(false);
                    ((ImageButton)findViewById(R.id.ib_settings)).setEnabled(false);
                    enbaleSettings(true);
                    break;
                }
                case R.id.ib_stop:{
                    Settings.setPlayback(true);

                    StopRecord();
                    break;
                }
                case R.id.b_save:{
                    enableButtons(false);
                    enbaleSettings(false);
                    updateSettings();
                    break;
                }
                case R.id.b_reset:{
                    resetSettings();
                    break;
                }
                case R.id.b_enchanced:{
                    if(isEnchanced){
                        ((Button) findViewById(R.id.b_enchanced)).setBackgroundColor(Color.RED);
                        isEnchanced=false;
                    }else {
                        ((Button) findViewById(R.id.b_enchanced)).setBackgroundColor(Color.GREEN);
                        isEnchanced=true;
                    }
                    enchancedDOA();
                    break;
                }
            }
        }
    };

    private void StopRecord() {
        switch(mode){
            case 1:{
                waveRecorder.stopRecording();
                waveRecorder = null;
                enableButtons(false);
                signalProcessing.release();
                break;
            }
        }
    }

    protected void resetSettings() {
        etFrameTime.setText(R.string.etFrametime);
        Settings.setupFrameTime(Settings.frameTime);
        etAngle.setText(R.string.etStepsize);
        Settings.setupStepSize(Settings.stepSize);
        Settings.setupThreadTime(Settings.ThreadTime);

        Settings.setupPlayback(1);
        Settings.setupNoiseType(1);
    }

    protected void updateSettings() {
        float initial_angle = 0;
        float initial_frame = 0;

        //Settings.setupStepSize(Float.parseFloat(etAngle.getText().toString()));
        Log.e("MainGUI","StepSize:"+Float.toString(Settings.stepFactor));
        //Settings.setupFrameTime(Float.parseFloat(etFrameTime.getText().toString()));
        Settings.setupThreadTime(Float.parseFloat(etFrameTime.getText().toString()));
        Settings.setupDurationTime(Float.parseFloat(etThreadDuration.getText().toString()));
        Log.e("MainGUI","FrameTime:"+Float.toString(Settings.frameTime));
        //myradiobutton = (RadioButton)findViewById(selectId);
        Settings.setupPlayback(playback);
        Log.e("MainGUI","Playback:"+Integer.toString(Settings.debugLevel));
    }

    private void StartRecord() {
        if (ActivityCompat.checkSelfPermission(this, Manifest.permission.RECORD_AUDIO) !=
                PackageManager.PERMISSION_GRANTED) {
            //   statusView.setText(getString(R.string.status_record_perm));
            ActivityCompat.requestPermissions(
                    this,
                    new String[] { Manifest.permission.RECORD_AUDIO },
                    AUDIO_ECHO_REQUEST);
            return;
        }

        if (ActivityCompat.checkSelfPermission(this, Manifest.permission.WRITE_EXTERNAL_STORAGE) !=
                PackageManager.PERMISSION_GRANTED) {
            //   statusView.setText(getString(R.string.status_record_perm));
            ActivityCompat.requestPermissions(
                    this,
                    new String[] { Manifest.permission.WRITE_EXTERNAL_STORAGE },
                    AUDIO_ECHO_REQUEST);
            return;
        }

        if (ActivityCompat.checkSelfPermission(this, Manifest.permission.MODIFY_AUDIO_SETTINGS) !=
                PackageManager.PERMISSION_GRANTED) {
            //   statusView.setText(getString(R.string.status_record_perm));
            ActivityCompat.requestPermissions(
                    this,
                    new String[] { Manifest.permission.MODIFY_AUDIO_SETTINGS },
                    AUDIO_ECHO_REQUEST);
            return;
        }

        inputFrames.clear();
        //inputFrames2.clear();
        processedFrames.clear();
//inputFrames.size();
        new Processing(inputFrames,processedFrames);


        waveRecorder = new WaveRecorder(inputFrames);

        if(Settings.debugLevel!=2){
            mode = 1;
            //new AudioSaver(getString(R.string.appDirectory)+ File.separator+System.currentTimeMillis(),processedFrames);
        }


        processingFlag.set(true);
        enableButtons(true);
    }
    private void enchancedDOA() {


        Settings.setupEnchanced(isEnchanced);

    }

    public int getMode() {
        return mode;
    }

    public void done() {
        if(processingFlag.getAndSet(false)) {
            runOnUiThread(
                    new Runnable() {
                        public void run() {
                            enableButtons(false);
                        }
                    }
            );
        }
    }

    public void notify(String message) {

    }
    /*
     * Loading our Libs
     */
    static {
        System.loadLibrary("echo");
    }



    public void onResume() {
        super.onResume();




        // setContentView(animation);



    }

    public void updateGUI()
    {
        //CircleView circle=new CircleView(this,45) ;
        circleanimation animation = new circleanimation(Settings.et_circle, Settings.output_angle);
        //animation.setInterpolator(new BounceInterpolator());
        //animation.setDuration(1000);
        Settings.et_circle.startAnimation(animation);

        String output_a = Float.toString(Settings.output_angle);
        
        Settings.et_degree.setText(output_a);

        String output_ab = Float.toString(Settings.SFx_avg);
        //String output_ab1 = output_ab.substring(0,3);
        //output_ab = String.format("%.4f", output_ab);
        Settings.et_SFxavg.setText(output_ab);
        Settings.ake_orientation.setTextSize(22);
        Settings.ake_orientation.setText("Re-orient by "+output_a+" degree ");
    }

    public void updateGUI_text()
    {
        float output_degree = 0;
        long tmp = 15;
        if(Settings.angle_counter<2*tmp)
        {
            output_degree = 0;
            Settings.angle_counter++;
        }
        else if(Settings.angle_counter<3*tmp+3)
        {
            output_degree = 15;
            Settings.angle_counter++;
        }
        else if(Settings.angle_counter<4*tmp)
        {
            output_degree = 30;
            Settings.angle_counter++;
        }
        else if(Settings.angle_counter<5*tmp-2)
        {
            output_degree = 45;
            Settings.angle_counter++;
        }
        else if(Settings.angle_counter<6*tmp+3)
        {
            output_degree = 60;
            Settings.angle_counter++;
        }
        else if(Settings.angle_counter<7*tmp)
        {
            output_degree = 75;
            Settings.angle_counter++;
        }
        else if(Settings.angle_counter<8*tmp-2)
        {
            output_degree = 90;
            Settings.angle_counter++;
        }
        else if(Settings.angle_counter<9*tmp)
        {
            output_degree = 105;
            Settings.angle_counter++;
        }
        else if(Settings.angle_counter<10*tmp)
        {
            output_degree = 120;
            Settings.angle_counter++;
        }
        else if(Settings.angle_counter<11*tmp)
        {
            output_degree = 135;
            Settings.angle_counter++;
        }
        else if(Settings.angle_counter<12*tmp)
        {
            output_degree = 150;
            Settings.angle_counter++;
        }
        else if(Settings.angle_counter<13*tmp)
        {
            output_degree = 165;
            Settings.angle_counter++;
        }
        else if(Settings.angle_counter<14*tmp)
        {
            output_degree = 180;
            Settings.angle_counter++;
        }
        else if(Settings.angle_counter<15*tmp-2)
        {
            output_degree = 165;
            Settings.angle_counter++;
        }
        else if(Settings.angle_counter<16*tmp)
        {
            output_degree = 150;
            Settings.angle_counter++;
        }
        else if(Settings.angle_counter<17*tmp)
        {
            output_degree = 135;
            Settings.angle_counter++;
        }
        else if(Settings.angle_counter<18*tmp)
        {
            output_degree = 120;
            Settings.angle_counter++;
        }
        else if(Settings.angle_counter<19*tmp)
        {
            output_degree = 105;
            Settings.angle_counter++;
        }
        else if(Settings.angle_counter<20*tmp)
        {
            output_degree = 90;
            Settings.angle_counter++;
        }
        else if(Settings.angle_counter<21*tmp+6)
        {
            output_degree = 75;
            Settings.angle_counter++;
        }
        else
        {
            Settings.angle_counter = 14*tmp;
        }

        Settings.angle_counter++;

        //CircleView circle=new CircleView(this,45) ;
        circleanimation animation = new circleanimation(Settings.et_circle, output_degree);

        Settings.et_circle.startAnimation(animation);

        String output_a = Float.toString(output_degree);

        Settings.et_degree.setText(output_a);

        //output_a = Float.toString(Settings.SFx_avg);
        Settings.et_SFxavg.setText(output_a);
    }



}
