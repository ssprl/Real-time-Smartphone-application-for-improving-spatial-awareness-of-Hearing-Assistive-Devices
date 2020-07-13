package com.example.axk166230.TwoChanDOA_v4;

/**
 * Created by axk166230 on 9/11/2017.
 */

public class SignalProcessing {

    private long pointer;
    public SignalProcessing(){
        pointer = paramInitilization(Settings.Fs,Settings.stepSize,Settings.frameSize, Settings.stepFactor, Settings.noisetype, Settings.ThreadTime, Settings.DurationTime,Settings.isEnchanced);
    }

    public void frameProcess(short[] input){
        realtimeProcessing(pointer,input);
    }

    public float getTime(){
        return getTime(pointer);
    }

    public float getComputeTime(){
        return getComputeTime(pointer);
    }

    public float getFilteringTime(){
        return getFilteringTime(pointer);
    }

    public void release(){
        paramElimination(pointer);
    }

    public short[] getSoundOutput(){
        return soundOutput(pointer,Settings.debugLevel);
    }

    public float[] getAngleOutput(){
        return angleOutput(pointer,Settings.debugLevel);
    }

    public short[] getDataOutput(){
        return dataOutput(pointer,Settings.debugLevel);
    }


    public static native long paramInitilization(int frequency, int stepSize,int frameSize,float stepFactor, int noisetype, float ThreadTime, float DurationTime, boolean isEnchanced);
    public static native void paramElimination(long memoryPointer);
    public static native float getTime(long memoryPointer);
    public static native float getComputeTime(long memoryPointer);
    public static native float getFilteringTime(long memoryPointer);
    public static native void realtimeProcessing(long memoryPointer,short[] input);
    public static native short[] soundOutput(long memoryPointer,int outputSelection);

    public static native float[] angleOutput(long memoryPointer,int outputSelection);

    public static native short[] dataOutput(long memoryPointer,int outputSelection);


}

