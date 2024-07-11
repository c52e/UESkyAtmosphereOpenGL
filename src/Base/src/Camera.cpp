#include "Camera.h"

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Utils.h"

Camera::Camera(glm::vec3 position, float yaw, float pitch, float aspect)
	: position_(position), up_(kWorldUp), yaw_(yaw), pitch_(pitch), aspect_(aspect) {
	UpdateVectors();
}

glm::mat4 Camera::ViewMatrix() const {
	return glm::lookAt(position_, position_ + front_, up_);
}

glm::mat4 Camera::ProjectionInvertedInfinity() const {
	float tanHalfFovy = glm::tan(glm::radians(fovy) * 0.5f);
	return glm::mat4{
		{1.0f / (aspect_ * tanHalfFovy), 0.0f, 0.0f, 0.0f},
		{0.0f, 1.0f / tanHalfFovy, 0.0f, 0.0f},
		{0.0f, 0.0f, 0.0f, -1.0f},
		{0.0f, 0.0f, zNear, 0.0f}
	};
}

glm::mat4 Camera::ProjectionFarClip(float zFar) const {
	return glm::perspective(glm::radians(fovy), aspect_, zNear, zFar);
}

void Camera::Rotate(float dPitch, float dYaw) {
	pitch_ += dPitch;
	yaw_ += dYaw;

	if (pitch_ > 89.0f)
		pitch_ = 89.0f;
	else if (pitch_ < -89.0f)
		pitch_ = -89.0f;

	UpdateVectors();
}

void Camera::UpdateVectors() {
	FromThetaPhiToDirection(glm::radians(90.f - pitch_), glm::radians(yaw_), glm::value_ptr(front_));
	right_ = glm::normalize(glm::cross(front_, kWorldUp));
	up_ = glm::normalize(glm::cross(right_, front_));
}
