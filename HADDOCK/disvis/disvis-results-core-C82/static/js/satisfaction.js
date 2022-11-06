
var ratingScore;

$(".rating-icon").click(function () {
    ratingScore = this.id.substring(6)
    // Hide rating icons
    $('.rating-icon').hide()
    // Show reCAPTCHA div
    $('.g-recaptcha').show()
});


function recaptchaCallback() {
    // grecaptcha object exists after this function is called
    function isCaptchaChecked() {
        return grecaptcha && grecaptcha.getResponse().length !== 0;
    }
    if (isCaptchaChecked()) {
        // Captcha has been clicked and validated, post!
        $.ajax({
            url: "/satisfaction",
            // url: "{{ request.host_url }}/satisfaction",
            type: 'POST',
            data: {
                rating: ratingScore,
                g_recaptcha_response: grecaptcha.getResponse()
            },
            async: true,
        }).done(function (data) {
            // Hide the reCAPTCHA div
            $('.g-recaptcha').hide()
            // Show the rated div
            $('#rated').show()
        });
    }
}
