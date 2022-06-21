//require jquery
//require util.js

Mouse3D = makeClass({
	/*
	args:
		pressObj: default window
		releaseObj: default window 
		mouseup
		mousedown
		move : function(dx,dy)  for click-and-drag  which match mouse.deltaX, mouse.deltaY
		passiveMove : function(dx,dy)  for simple movement
		zoom : function(dz)
		click : function(event)
		touchclickstart
		touchclickend
		preventDefault - explicit set to false or all defaults will be prevented

	params:
		lastX
		lastY
		downX
		downY
		xf
		yf
		isDown
		isDragging
		isTouchDown
		zoomTouchPts
		newZoomTouchPts

	*/
	init : function(args) {
		this.xf = .5;
		this.yf = .5;
		this.isDown = false;
		this.isDragging = false;
		this.isTouchDown = false;
		this.zoomTouchPts = [[0,0],[0,0]];
		this.newZoomTouchPts = [[0,0],[0,0]];
		this.preventDefault = !(args.preventDefault === false);
		
		var thiz = this;
		this.pressObj = args.pressObj !== undefined ? args.pressObj : window;
		this.releaseObj = args.releaseObj !== undefined ? args.releaseObj : window;
		var callbackNames = ['mouseup','mousedown','move','passiveMove','zoom','click','touchclickstart','touchclickend'];
		$.each(callbackNames, function(_,callbackName) {
			if (args[callbackName] !== undefined) {
				thiz[callbackName] = args[callbackName];
			}
		});

		//keep this member around for unbinding in the event of finding touch events are used
		//therefore don't expect 'this' to be referenced correctly
		this.mousedownCallback = function(e) {
			if (thiz.preventDefault) e.preventDefault();
			thiz.doMouseDown(e);
		};
		$(this.pressObj).bind('mousedown', this.mousedownCallback);
	
		this.mouseupCallback = function(e) {
			if (thiz.preventDefault) e.preventDefault();
			thiz.doMouseUp(e);
		};
		$(this.releaseObj).bind('mouseup', this.mouseupCallback);
		
		this.mousemoveCallback = function(e) {
			if (thiz.preventDefault) e.preventDefault();
			thiz.doMouseMove(e);
		};
		$(this.pressObj).bind('mousemove', this.mousemoveCallback);
		
		this.mousewheelCallback = function(e) {
			if (thiz.preventDeafult) e.preventDefault();
			var zoomChange = e.originalEvent.wheelDelta;
			if (thiz.zoom) {
				thiz.zoom(zoomChange, 'wheel');
				e.preventDefault();
			}
		};
		$(this.pressObj).bind('mousewheel', this.mousewheelCallback);
		
		//special case for Firefox 25:
		//http://www.javascriptkit.com/javatutors/onmousewheel.shtml
		$(this.pressObj).bind('DOMMouseScroll', function(e) {
			if (thiz.preventDefault) e.preventDefault();
			var zoomChange = e.originalEvent.detail;
			if (thiz.zoom) thiz.zoom(zoomChange * -120, 'wheel');
		});
		
		this.clickCallback = function(e) {
			//TODO also check l-infinite distance?  or total distance while mousedown travelled?
			if (thiz.click) thiz.click(e);
		};
		$(this.pressObj).bind('click', this.clickCallback);

		var unbindMouse = function() {
			$(thiz.pressObj).unbind('mousedown', thiz.mousedownCallback);
			$(thiz.releaseObj).unbind('mouseup', thiz.mouseupCallback);
			$(thiz.pressObj).unbind('mousemove', thiz.mousemoveCallback);
			$(thiz.pressObj).unbind('mousewheel', thiz.mousewheelCallback);
			$(thiz.pressObj).unbind('click', thiz.clickCallback);
		};

		$(this.pressObj)
			.bind('touchstart', function(e) { 	
				unbindMouse();
				thiz.isTouchDown = false;
				if (thiz.preventDefault) e.preventDefault();
				thiz.doMouseDown(e.originalEvent.targetTouches[0]);
				if (thiz.touchclickstart) thiz.touchclickstart();
			})
			.bind('touchmove', function(e) {
				if (thiz.preventDefault) e.preventDefault();
				if (e.originalEvent.touches.length >= 2) {
					//do a pinch zoom
					if (!thiz.isTouchDown) {
						//record current events
						thiz.getTouchPts(e, thiz.zoomTouchPts);
						thiz.isTouchDown = true;
					} else {
						//do zoom based on movement
						thiz.getTouchPts(e, thiz.newZoomTouchPts);
						var zoomChange = thiz.calcDist(thiz.newZoomTouchPts) - thiz.calcDist(thiz.zoomTouchPts);
						if (zoomChange != 0) {
							if (thiz.zoom) {
								thiz.zoom(100 * zoomChange, 'pinch');
								thiz.zoomTouchPts[0][0] = thiz.newZoomTouchPts[0][0];
								thiz.zoomTouchPts[0][1] = thiz.newZoomTouchPts[0][1];
								thiz.zoomTouchPts[1][0] = thiz.newZoomTouchPts[1][0];
								thiz.zoomTouchPts[1][1] = thiz.newZoomTouchPts[1][1];
							}
						}
					}
				} else {
					//don't reset zoom once we've begun touch zooming
					//until we're finished with the gesture
					//thiz.isTouchDown = false;
				}
				if (!thiz.isTouchDown) {	//only rotate if we haven't begun zooming
					thiz.doMouseMove(e.originalEvent.targetTouches[0]);
				}
			})
			.bind('touchend touchcancel', function(e) {
				if (thiz.preventDefault) e.preventDefault();
				var touch = e.originalEvent.changedTouches[0];
				var upPosX = touch.pageX;
				var upPosY = touch.pageY;
				thiz.deltaX = upPosX - thiz.downX;
				thiz.deltaY = upPosY - thiz.downY;
				thiz.xf = upPosX / window.innerWidth;
				thiz.yf = upPosY / window.innerHeight;
				var linf = Math.max( Math.abs(thiz.deltaX), Math.abs(thiz.deltaY) );
				if (linf < 2) {
					if (thiz.touchclickend) thiz.touchclickend();
					if (thiz.click) thiz.click(touch);
				} 
				thiz.doMouseUp(touch);
			})
		;
	},
	doMouseDown : function(e) {
		this.isDragging = false;
		this.isDown = true;
		this.lastX = e.pageX;
		this.lastY = e.pageY;
		this.xf = e.pageX / window.innerWidth;
		this.yf = e.pageY / window.innerHeight;
		this.downX = this.lastX;
		this.downY = this.lastY;
		this.deltaX = 0;
		this.deltaY = 0;

		if (this.mousedown) this.mousedown(e);
	},
	doMouseUp : function(e) {
		this.isDown = false;
		if (this.mouseup) this.mouseup();
	},
	doMouseMove : function(e) {
		var thisX = e.pageX;
		var thisY = e.pageY;
		this.xf = thisX / window.innerWidth;
		this.yf = thisY / window.innerHeight;

		this.deltaX = thisX - this.lastX;
		this.deltaY = thisY - this.lastY;
		this.lastX = thisX;
		this.lastY = thisY;

		if (this.isDown) {
			this.isDragging = true;
			if (e.shiftKey) {
				if (this.zoom) this.zoom(-100 * this.deltaY, 'shift');
			} else {
				if (this.move) this.move(this.deltaX, this.deltaY);
			}
		} else {
			if (this.passiveMove) this.passiveMove(this.deltaX, this.deltaY);
		}
	},
	getTouchPts : function(e, pts) {
		pts[0][0] = e.originalEvent.changedTouches[0].pageX;
		pts[0][1] = e.originalEvent.changedTouches[0].pageY;
		pts[1][0] = e.originalEvent.changedTouches[1].pageX;
		pts[1][1] = e.originalEvent.changedTouches[1].pageY;
	},
	calcDist : function(pts) {
		var dx = pts[0][0] - pts[1][0];
		var dy = pts[0][1] - pts[1][1];
		return Math.sqrt(dx*dx + dy*dy);
	}
});

