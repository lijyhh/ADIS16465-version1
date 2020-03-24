/*
	
ap_ahrs_quaternion代码，基于杰布马德威克四元数码

看到http://www.x-io.co.uk/res/doc/madgwick_internal_report.pdf

适于APM安得烈Tridgell基于最初的想法，

从贾斯廷的山毛榉的讨论和原型。

这个库是自由软件；你可以重新分配和/或

修改它的GNU通用公共条件下

许可证由自由软件基金会发布的；或

2.1版的许可证，或（根据你的选择）任何时候

版。
*/	
#include <FastSerial.h>
#include <AP_AHRS.h>

// to keep the code as close to the original as possible, we use these
// macros for quaternion access
#define SEq_1 q.q1
#define SEq_2 q.q2
#define SEq_3 q.q3
#define SEq_4 q.q4

// 函数来计算一个数的迭代没有磁强计

void AP_AHRS_Quaternion::update_IMU(float deltat, Vector3f &gyro, Vector3f &accel)
{
    // Local system variables
    float norm;                                                            // vector norm
    float SEqDot_omega_1, SEqDot_omega_2, SEqDot_omega_3, SEqDot_omega_4; // quaternion derrivative from gyroscopes elements
    float f_1, f_2, f_3;                                                   // objective function elements
    float J_11or24, J_12or23, J_13or22, J_14or21, J_32, J_33;              // objective function Jacobian elements
    float SEqHatDot_1, SEqHatDot_2, SEqHatDot_3, SEqHatDot_4;              // estimated direction of the gyroscope error

    // Axulirary variables to avoid reapeated calcualtions
    float halfSEq_1 = 0.5f * SEq_1;
    float halfSEq_2 = 0.5f * SEq_2;
    float halfSEq_3 = 0.5f * SEq_3;
    float halfSEq_4 = 0.5f * SEq_4;
    float twoSEq_1 = 2.0f * SEq_1;
    float twoSEq_2 = 2.0f * SEq_2;
    float twoSEq_3 = 2.0f * SEq_3;

    // 	估计的陀螺误差方向（弧度）
    
    Vector3f w_err;

    // 	正常化加速度矢量
    
    accel.normalize();
    if (accel.is_inf()) {
	    // discard this data point
	    renorm_range_count++;
	    return;
    }

    // 	计算目标函数和雅可比矩阵
    
    f_1 = twoSEq_2 * SEq_4 - twoSEq_1 * SEq_3 - accel.x;
    f_2 = twoSEq_1 * SEq_2 + twoSEq_3 * SEq_4 - accel.y;
    f_3 = 1.0f - twoSEq_2 * SEq_2 - twoSEq_3 * SEq_3 - accel.z;
    J_11or24 = twoSEq_3;                                     // J_11 negated in matrix multiplication
    J_12or23 = 2.0f * SEq_4;
    J_13or22 = twoSEq_1;                                     // J_12 negated in matrix multiplication
    J_14or21 = twoSEq_2;
    J_32 = 2.0f * J_14or21;                                  // negated in matrix multiplication
    J_33 = 2.0f * J_11or24;                                  // negated in matrix multiplication

    // 	梯度的计算（矩阵乘法）
    
    SEqHatDot_1 = J_14or21 * f_2 - J_11or24 * f_1;
    SEqHatDot_2 = J_12or23 * f_1 + J_13or22 * f_2 - J_32 * f_3;
    SEqHatDot_3 = J_12or23 * f_2 - J_33 * f_3 - J_13or22 * f_1;
    SEqHatDot_4 = J_14or21 * f_1 + J_11or24 * f_2;

    // 	正常化的梯度
    
    norm = 1.0/safe_sqrt(SEqHatDot_1 * SEqHatDot_1 + SEqHatDot_2 * SEqHatDot_2 + SEqHatDot_3 * SEqHatDot_3 + SEqHatDot_4 * SEqHatDot_4);
    if (isinf(norm)) {
	    // we can't do an update - discard this data point and
	    // hope the next one is better
	    renorm_range_count++;
	    return;
    }
    SEqHatDot_1 *= norm;
    SEqHatDot_2 *= norm;
    SEqHatDot_3 *= norm;
    SEqHatDot_4 *= norm;

    // 	计算四元数derrivative的测量陀螺仪
    
    SEqDot_omega_1   = -halfSEq_2 * gyro.x - halfSEq_3 * gyro.y - halfSEq_4 * gyro.z;
    SEqDot_omega_2   = halfSEq_1  * gyro.x + halfSEq_3 * gyro.z - halfSEq_4 * gyro.y;
    SEqDot_omega_3   = halfSEq_1  * gyro.y - halfSEq_2 * gyro.z + halfSEq_4 * gyro.x;
    SEqDot_omega_4   = halfSEq_1  * gyro.z + halfSEq_2 * gyro.y - halfSEq_3 * gyro.x;

    //	计算然后将估计的四元数derrivative
    
    SEq_1 += (SEqDot_omega_1 - (beta * SEqHatDot_1)) * deltat;
    SEq_2 += (SEqDot_omega_2 - (beta * SEqHatDot_2)) * deltat;
    SEq_3 += (SEqDot_omega_3 - (beta * SEqHatDot_3)) * deltat;
    SEq_4 += (SEqDot_omega_4 - (beta * SEqHatDot_4)) * deltat;

    // Normalise quaternion
    norm = 1.0/safe_sqrt(SEq_1 * SEq_1 + SEq_2 * SEq_2 + SEq_3 * SEq_3 + SEq_4 * SEq_4);
    if (isinf(norm)) {
	    // our quaternion is bad! Reset based on roll/pitch/yaw
	    // and hope for the best ...
	    renorm_blowup_count++;
	    if (_compass) {
		    _compass->null_offsets_disable();
	    }
	    q.from_euler(roll, pitch, yaw);
	    if (_compass) {
		    _compass->null_offsets_enable();
	    }
	    return;
    }
    SEq_1 *= norm;
    SEq_2 *= norm;
    SEq_3 *= norm;
    SEq_4 *= norm;
}


// 函数来计算一个数的迭代包括磁强计

void AP_AHRS_Quaternion::update_MARG(float deltat, Vector3f &gyro, Vector3f &accel, Vector3f &mag)
{
    // local system variables
    float norm;                                                             //                vector norm
    float SEqDot_omega_1, SEqDot_omega_2, SEqDot_omega_3, SEqDot_omega_4;   //             从陀螺仪的元素数率
    
    float f_1, f_2, f_3, f_4, f_5, f_6;                                     //               目标函数的元素
    
    float J_11or24, J_12or23, J_13or22, J_14or21, J_32, J_33,               //                目标函数的雅可比矩阵元素
    
    J_41, J_42, J_43, J_44, J_51, J_52, J_53, J_54, J_61, J_62, J_63, J_64; //
    float SEqHatDot_1, SEqHatDot_2, SEqHatDot_3, SEqHatDot_4;               //              	估计的陀螺误差方向
    

    // 	计算的通量在地球的框架
    
    Vector3f flux;

    // 	估计的陀螺误差方向（弧度）
    
    Vector3f w_err;

    // 	正常化加速度矢量
    
    accel.normalize();
    if (accel.is_inf()) {
	    // 		丢弃这个数据点
	    
	    renorm_range_count++;
	    return;
    }

    // 	正常化的磁强计测量
    
    mag.normalize();
    if (mag.is_inf()) {
	    //丢弃这个数据点
	    
	    renorm_range_count++;
	    return;
    }

    //	为了避免重复计算的辅助变量
    
    float halfSEq_1 = 0.5f * SEq_1;
    float halfSEq_2 = 0.5f * SEq_2;
    float halfSEq_3 = 0.5f * SEq_3;
    float halfSEq_4 = 0.5f * SEq_4;
    float twoSEq_1 = 2.0f * SEq_1;
    float twoSEq_2 = 2.0f * SEq_2;
    float twoSEq_3 = 2.0f * SEq_3;
    float twoSEq_4 = 2.0f * SEq_4;
    float twob_x = 2.0f * b_x;
    float twob_z = 2.0f * b_z;
    float twob_xSEq_1 = 2.0f * b_x * SEq_1;
    float twob_xSEq_2 = 2.0f * b_x * SEq_2;
    float twob_xSEq_3 = 2.0f * b_x * SEq_3;
    float twob_xSEq_4 = 2.0f * b_x * SEq_4;
    float twob_zSEq_1 = 2.0f * b_z * SEq_1;
    float twob_zSEq_2 = 2.0f * b_z * SEq_2;
    float twob_zSEq_3 = 2.0f * b_z * SEq_3;
    float twob_zSEq_4 = 2.0f * b_z * SEq_4;
    float SEq_1SEq_2;
    float SEq_1SEq_3 = SEq_1 * SEq_3;
    float SEq_1SEq_4;
    float SEq_2SEq_3;
    float SEq_2SEq_4 = SEq_2 * SEq_4;
    float SEq_3SEq_4;
    Vector3f twom = mag * 2.0;

    // 	计算目标函数和雅可比矩阵
    
    f_1 = twoSEq_2 * SEq_4 - twoSEq_1 * SEq_3 - accel.x;
    f_2 = twoSEq_1 * SEq_2 + twoSEq_3 * SEq_4 - accel.y;
    f_3 = 1.0f - twoSEq_2 * SEq_2 - twoSEq_3 * SEq_3 - accel.z;
    f_4 = twob_x * (0.5f - SEq_3 * SEq_3 - SEq_4 * SEq_4) + twob_z * (SEq_2SEq_4 - SEq_1SEq_3) - mag.x;
    f_5 = twob_x * (SEq_2 * SEq_3 - SEq_1 * SEq_4) + twob_z * (SEq_1 * SEq_2 + SEq_3 * SEq_4) - mag.y;
    f_6 = twob_x * (SEq_1SEq_3 + SEq_2SEq_4) + twob_z * (0.5f - SEq_2 * SEq_2 - SEq_3 * SEq_3) - mag.z;
    J_11or24 = twoSEq_3;                                  // 11否定矩阵乘法
    J_12or23 = 2.0f * SEq_4;
    J_13or22 = twoSEq_1;                                  // J_12 negated in matrix multiplication
    J_14or21 = twoSEq_2;
    J_32 = 2.0f * J_14or21;                               // 	否定的矩阵乘法
    
    J_33 = 2.0f * J_11or24;                               // negated in matrix multiplication
    J_41 = twob_zSEq_3;                                   // negated in matrix multiplication
    J_42 = twob_zSEq_4;
    J_43 = 2.0f * twob_xSEq_3 + twob_zSEq_1;              // negated in matrix multiplication
    J_44 = 2.0f * twob_xSEq_4 - twob_zSEq_2;              // negated in matrix multiplication
    J_51 = twob_xSEq_4 - twob_zSEq_2;                     // negated in matrix multiplication
    J_52 = twob_xSEq_3 + twob_zSEq_1;
    J_53 = twob_xSEq_2 + twob_zSEq_4;
    J_54 = twob_xSEq_1 - twob_zSEq_3;                     // negated in matrix multiplication
    J_61 = twob_xSEq_3;
    J_62 = twob_xSEq_4 - 2.0f * twob_zSEq_2;
    J_63 = twob_xSEq_1 - 2.0f * twob_zSEq_3;
    J_64 = twob_xSEq_2;

    // 	梯度的计算（矩阵乘法）
    
    SEqHatDot_1 = J_14or21 * f_2 - J_11or24 * f_1              - J_41 * f_4 - J_51 * f_5 + J_61 * f_6;
    SEqHatDot_2 = J_12or23 * f_1 + J_13or22 * f_2 - J_32 * f_3 + J_42 * f_4 + J_52 * f_5 + J_62 * f_6;
    SEqHatDot_3 = J_12or23 * f_2 - J_33 * f_3 - J_13or22 * f_1 - J_43 * f_4 + J_53 * f_5 + J_63 * f_6;
    SEqHatDot_4 = J_14or21 * f_1 + J_11or24 * f_2              - J_44 * f_4 - J_54 * f_5 + J_64 * f_6;

    // 	正常化梯度对陀螺误差方向估计
    
    norm = 1.0 / safe_sqrt(SEqHatDot_1 * SEqHatDot_1 + SEqHatDot_2 * SEqHatDot_2 + SEqHatDot_3 * SEqHatDot_3 + SEqHatDot_4 * SEqHatDot_4);
    if (isinf(norm)) {
	    // discard this data point
	    renorm_range_count++;
	    return;
    }
    SEqHatDot_1 *= norm;
    SEqHatDot_2 *= norm;
    SEqHatDot_3 *= norm;
    SEqHatDot_4 *= norm;

    // 	计算陀螺仪误差角估计的方向
    
    w_err.x = twoSEq_1 * SEqHatDot_2 - twoSEq_2 * SEqHatDot_1 - twoSEq_3 * SEqHatDot_4 + twoSEq_4 * SEqHatDot_3;
    w_err.y = twoSEq_1 * SEqHatDot_3 + twoSEq_2 * SEqHatDot_4 - twoSEq_3 * SEqHatDot_1 - twoSEq_4 * SEqHatDot_2;
    w_err.z = twoSEq_1 * SEqHatDot_4 - twoSEq_2 * SEqHatDot_3 + twoSEq_3 * SEqHatDot_2 - twoSEq_4 * SEqHatDot_1;

    //	跟踪误差率
    
    _error_rp_sum += 0.5*(fabs(w_err.x) + fabs(w_err.y));
    _error_yaw_sum += fabs(w_err.z);
    _error_rp_count++;
    _error_yaw_count++;

    // 	计算陀螺仪偏置三角洲
    
    Vector3f drift_delta = w_err * (deltat * zeta);

  
  //不允许漂移率将超过。这防止了 
  //突然漂移变化的罗盘来自一个中断
  
    float max_change = _gyro_drift_limit * deltat;
    drift_delta.x = constrain(drift_delta.x, -max_change, max_change);
    drift_delta.y = constrain(drift_delta.y, -max_change, max_change);
    drift_delta.z = constrain(drift_delta.z, -max_change, max_change);
    gyro_bias += drift_delta;

    // 	正确的陀螺漂移的阅读
    
    gyro -= gyro_bias;

    // 	计算四元数率的测量陀螺仪
    
    SEqDot_omega_1 = -halfSEq_2 * gyro.x - halfSEq_3 * gyro.y - halfSEq_4 * gyro.z;
    SEqDot_omega_2 = halfSEq_1  * gyro.x + halfSEq_3 * gyro.z - halfSEq_4 * gyro.y;
    SEqDot_omega_3 = halfSEq_1  * gyro.y - halfSEq_2 * gyro.z + halfSEq_4 * gyro.x;
    SEqDot_omega_4 = halfSEq_1  * gyro.z + halfSEq_2 * gyro.y - halfSEq_3 * gyro.x;

    // 	然后将估计的四元数率的计算
    
    SEq_1 += (SEqDot_omega_1 - (beta * SEqHatDot_1)) * deltat;
    SEq_2 += (SEqDot_omega_2 - (beta * SEqHatDot_2)) * deltat;
    SEq_3 += (SEqDot_omega_3 - (beta * SEqHatDot_3)) * deltat;
    SEq_4 += (SEqDot_omega_4 - (beta * SEqHatDot_4)) * deltat;

    // 	正常化的四元数
    
    norm = 1.0/safe_sqrt(SEq_1 * SEq_1 + SEq_2 * SEq_2 + SEq_3 * SEq_3 + SEq_4 * SEq_4);
    if (isinf(norm)) {
	    
		//我们的四元数差！复位基于滚/俯仰/偏航
		//和最好的希望…
	    
	    renorm_blowup_count++;
	    _compass->null_offsets_disable();
	    q.from_euler(roll, pitch, yaw);
	    _compass->null_offsets_disable();
	    return;
    }
    SEq_1 *= norm;
    SEq_2 *= norm;
    SEq_3 *= norm;
    SEq_4 *= norm;

    
	//计算通量在地球的框架
	//重新计算axulirary变量
    
    SEq_1SEq_2 = SEq_1 * SEq_2;
    SEq_1SEq_3 = SEq_1 * SEq_3;
    SEq_1SEq_4 = SEq_1 * SEq_4;
    SEq_3SEq_4 = SEq_3 * SEq_4;
    SEq_2SEq_3 = SEq_2 * SEq_3;
    SEq_2SEq_4 = SEq_2 * SEq_4;
    flux.x = twom.x * (0.5f - SEq_3 * SEq_3 - SEq_4 * SEq_4) + twom.y * (SEq_2SEq_3 - SEq_1SEq_4) + twom.z * (SEq_2SEq_4 + SEq_1SEq_3);
    flux.y = twom.x * (SEq_2SEq_3 + SEq_1SEq_4) + twom.y * (0.5f - SEq_2 * SEq_2 - SEq_4 * SEq_4) + twom.z * (SEq_3SEq_4 - SEq_1SEq_2);
    flux.z = twom.x * (SEq_2SEq_4 - SEq_1SEq_3) + twom.y * (SEq_3SEq_4 + SEq_1SEq_2) + twom.z * (0.5f - SEq_2 * SEq_2 - SEq_3 * SEq_3);

    // 	正常化的磁通矢量在X和Z的组件
    
    b_x = sqrt((flux.x * flux.x) + (flux.y * flux.y));
    b_z = flux.z;
}


// 函数来计算一个数的迭代

void AP_AHRS_Quaternion::update(void)
{
	Vector3f gyro, accel;
	float deltat;

	_imu->update();

	deltat = _imu->get_delta_time();
	if (deltat > 1.0) {
		
		//如果我们停止更新为1s，我们应该放弃这
		//输入数据。这可以在你运行的发生
		
		//代码调试器下，利用这个数据点
		
		//会抛弃你的态度由一个巨大的数额
		
		return;
	}

	if (!_have_initial_yaw && _compass &&
	    _compass->use_for_yaw()) {
		// setup the quaternion with initial compass yaw
		_compass->null_offsets_disable();
		q.from_euler(0, 0, _compass->heading);
		_have_initial_yaw = true;
		_compass_last_update = _compass->last_update;
		gyro_bias.zero();
		_compass->null_offsets_enable();
	}

	// get current IMU state
	gyro = _imu->get_gyro();

	//	四元数系统使用符号相反的加速
	
	accel = - _imu->get_accel();

	if (_centripetal && _gps && _gps->status() == GPS::GPS_OK) {
		// compensate for linear acceleration. This makes a
		// surprisingly large difference in the pitch estimate when
		// turning, plus on takeoff and landing
		float acceleration = _gps->acceleration();
		accel.x += acceleration;

		// 补偿向心加速度
		
		float veloc;
		veloc = _gps->ground_speed * 0.01;
		
		//小心本文计算的迹象。的
		
		//四元系统使用不同的迹象比
		
		//休息APM
		
		accel.y += (gyro.z - gyro_bias.z) * veloc;
		accel.z -= (gyro.y - gyro_bias.y) * veloc;
	}

	if (_compass != NULL && _compass->use_for_yaw()) {
		Vector3f mag = Vector3f(_compass->mag_x, _compass->mag_y, _compass->mag_z);
		update_MARG(deltat, gyro, accel, mag);
	} else {
		// 使用陀螺和accels四元数的解决步骤
		
		gyro -= gyro_bias;
		update_IMU(deltat, gyro, accel);
	}

#ifdef DESKTOP_BUILD
	if (q.is_nan()) {
		SITL_debug("QUAT NAN: deltat=%f roll=%f pitch=%f yaw=%f q=[%f %f %f %f] a=[%f %f %f] g=(%f %f %f)\n",
			   deltat, roll, pitch, yaw,
			   q.q1, q.q2, q.q3, q.q4,
			   accel.x, accel.y, accel.z,
			   gyro.x, gyro.y, gyro.z);
	}
#endif

	// 	保持校正陀螺的报告
	
	_gyro_corrected = gyro;

	// calculate our euler angles for high level control and navigation
	q.to_euler(&roll, &pitch, &yaw);

	// the code above assumes zero magnetic declination, so offset
	// the yaw here
	if (_compass != NULL) {
		yaw += _compass->get_declination();
	}

	// and integer Eulers
	roll_sensor  = 100 * ToDeg(roll);
	pitch_sensor = 100 * ToDeg(pitch);
	yaw_sensor   = 100 * ToDeg(yaw);
	if (yaw_sensor < 0) {
		yaw_sensor += 36000;
	}

}

// average error in roll/pitch since last call
float AP_AHRS_Quaternion::get_error_rp(void)
{
	float ret;
	if (_error_rp_count == 0) {
		return 0;
	}
	ret = _error_rp_sum / _error_rp_count;
	_error_rp_sum = 0;
	_error_rp_count = 0;
	return ret;
}

// average error in yaw since last call
float AP_AHRS_Quaternion::get_error_yaw(void)
{
	float ret;
	if (_error_yaw_count == 0) {
		return 0;
	}
	ret = _error_yaw_sum / _error_yaw_count;
	_error_yaw_sum = 0;
	_error_yaw_count = 0;
	return ret;
}

// reset attitude system
void AP_AHRS_Quaternion::reset(bool recover_eulers)
{
	if (recover_eulers) {
		q.from_euler(roll, pitch, yaw);
	} else {
		q(1, 0, 0, 0);
	}
	gyro_bias.zero();

	// reference direction of flux in earth frame
	b_x = 0;
	b_z = -1;
}
