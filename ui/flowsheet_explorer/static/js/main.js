function init()
{
    // Set up the event listeners for the slider that shows and hides the buttons
    var bubble = document.getElementById('bubble');
    var slider = document.getElementById('slider');

    slider.addEventListener('mousemove', moveBubble);
    slider.addEventListener('mousedown', show);
    slider.addEventListener('mouseup', hide);

    // var back_arrow = document.getElementById('back_arrow');
    // back_arrow.addEventListener('click', back_click);
}

$(function() {
    // This connects the back button with changing the images
    var form = $('form');
    $('#back_arrow').click(function(){
        var sliderVal = slider.value;
        if (sliderVal - 1 >= 0)
        {                
            $.ajax({
                type: "POST",
                url: form.action,
                data: $('#slider').val(parseInt($('#slider').val()) - 1),
                // Set the slider value to -= 1
                success: function(data){
                    $('#data').html(data);
                    init();
                    // Return the data html (all of it right now)
                    // Also call init so the slider bubble works
                }
            }).done(function(res){
            });
        }
    });
});

$(function() {
    // This connects the forward button with changing the images
    var form = $('form');
    $('#forward_arrow').click(function(){
        var sliderVal = slider.value;
        if ((parseInt(sliderVal) + 1) < range)
        {                
            $.ajax({
                type: "POST",
                url: form.action,
                data: $('#slider').val(parseInt($('#slider').val()) + 1),
                // Set the slider value to += 1
                success: function(data){
                    $('#data').html(data);
                    init();
                    // Return the data html (all of it right now)
                    // Also call init so the slider bubble works
                }
            }).done(function(res){
            });
        }
    });
});

$(function() {
    // This connects the slider changing to changing the images as well
    var form = $('form');
    $('#slider').on('change mouseup', function(){
        $.ajax({
            type: "POST",
            url: form.action,
            data: $(this).serialize(),
            success: function(data){
                $('#data').html(data);
                init();
                // Return the data html (all of it right now)
                // Also call init so the slider bubble works
            }
        }).done(function(res){
            //do something with the response from the server
        });
    });
});



var oldSliderVal = -1;   

var show = function(e)
{
    // Show the hover bubble
    bubble.style.left = e.clientX-(bubble.offsetWidth)+'px';  
    bubble.style.opacity = '1';
}

var hide = function()
{
    // Hide the hover bubble
    bubble.style.opacity = '0';
}

var moveBubble = function(e)
{
    // Move the bubble to follow the slider
    if(oldSliderVal !== '0' && oldSliderVal !== range)
    { 
        bubble.style.left = e.clientX-(bubble.offsetWidth)+'px';        
    }
    var sliderVal = slider.value
    bubble.innerHTML = slider_labels[sliderVal];
    oldSliderVal = sliderVal;
}
