package com.example.axk166230.TwoChanDOA_v4;

import android.view.animation.Animation;
import android.view.animation.Transformation;

/**
 * Created by axk166230 on 9/11/2017.
 */

public class circleanimation extends Animation {

    private Circle circle;

    private float oldAngle;
    private float newAngle;

    public circleanimation(Circle circle, float newAngle) {
        this.oldAngle = circle.getAngle();
        this.newAngle = newAngle;
        this.circle = circle;
    }

    @Override
    public void applyTransformation(float interpolatedTime, Transformation transformation) {
        float angle = oldAngle + ((newAngle - oldAngle) * interpolatedTime);

        //circle.;
        circle.setAngle(angle);
        circle.requestLayout();
    }




}


