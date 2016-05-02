/**
 * Trajectory control Widget for IPython
 *
 */
 require.config({
    "paths": {
      "jquery": "//ajax.googleapis.com/ajax/libs/jquery/2.0.0/jquery.min",
      "jqueryui": "nbextensions/jquery-ui.min"
    }
});

define(["widgets/js/widget",
        'jquery',
        'jqueryui'],
function( widget, $ ) {
    console.log('hello@@@');
    var TrajectoryControls = widget.DOMWidgetView.extend({

        render : function() {
            this.fps = this.model.get('fps');
            this.width = 600;
            this.height = 20;
            var model = this.model;
            // Create the ui elements
            //
            var tc = $('<div/>');
            tc.addClass("ui-widget-header")
              .addClass("ui-corner-all")
              .css({ padding: "6px", display: 'flex'});
            tc.width(this.width);

            var that = this;
            var slider = $('<div/>').slider({
                value: 0,
                max: this.model.get('n_frames'),
                slide: function (event, ui) {
                    model.set('frame', ui.value);
                    that.touch();
                },
                start: function(event, ui) {
                    that.pause();
                    that.running = false;
                },
            });
            this.running = false;

            var playButton = $('<button/>').button({
                    text: false,
                    icons: {
                        primary: "ui-icon-play"
                    },
                });

            playButton.click(function(event) {
                    if (!that.running) {
                        that.play();
                        that.running = true;
                    } else {
                        that.pause();
                        that.running = false;
                    }
                });
            playButton.width(this.height).height(this.height);
            playButton.css("float", "left");
            playButton.appendTo(tc);

            // Calculate the number of character to prepare for the correct space

            var maxLength = String(model.get("n_frames")).length;

            var frameIndicator = $("<span/>").text(model.get("frame") + "/" + model.get("n_frames"))
                                             .css( { "width": (2*maxLength) + "em",
                                                     "text-align": "right",
                                                     "margin-right": "6px" } );


            model.on("change:frame", function () {
                frameIndicator.text(model.get("frame") + "/" + model.get("n_frames"));
            });
            // We add an extra container for the slider just for the styling
            var sliderContainer = $("<div/>").css({ "margin": "4px 16px" ,
                                                    "flex-grow" : "1" });

            slider.appendTo(sliderContainer);
            sliderContainer.appendTo(tc);
            frameIndicator.appendTo(tc);

            // HACK: We need an extra layer because Jupyter widget creation is
            // stupid
            var dummy = $('<div/>');
            tc.appendTo(dummy);
            this.setElement(dummy);
            this.slider = slider;
            this.playButton = playButton;
        },


        update : function () {
            this.fps = this.model.get('fps');

            return TrajectoryControls.__super__.update.apply(this);
        },

        play : function () {
            var slider = this.slider;
            var that = this;
            this.playButton.button('option', { icons: { primary : "ui-icon-pause"} });

            if (slider.slider('value') == slider.slider('option', 'max')) {
                slider.slider('value', 0);
            }

            this.playCallbackId = setInterval( function () {

                if (slider.slider('value') < slider.slider('option', 'max')) {
                    slider.slider('value', slider.slider('value') + 1);
                    that.model.set('frame', slider.slider('value'));
                    that.touch();
                } else {
                    that.pause();
                    that.running = false;
                }

            }, 1000/this.fps);

        },

        pause : function () {
            this.playButton.button('option', { icons: { primary : "ui-icon-play"} });
            clearInterval(this.playCallbackId);
        },

        /* We receive custom messages from our python conterpart with DOMWidget.send */
        on_msg: function(msg) {
            if (msg.type == 'callMethod') {
                this[msg.methodName].call(this, msg.args);
            }

            return TrajectoryControls.__super__.update.apply(this);
        },

    });


    return {
        TrajectoryControls : TrajectoryControls
    };
});
