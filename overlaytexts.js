var showNames = true;

var overlayTexts = new function() {
	this.overlays = [];
	this.index = 0;
	this.maxOverlays = 20;

	this.updateBegin = function() {
		this.index = 0;
	};

	this.updateEnd = function() {
		for (var i = this.index; i < this.overlays.length; ++i) {
			var div = this.overlays[i].div
			var parent = div.parentNode;
			if (parent) parent.removeChild(div);
		}
	};

	this.add = function(target) {
		if (!showNames) return;
		
		for (var i = 0; i < this.index; ++i) {
			if (this.overlays[i].target == target) return;
		}
		
		var overlayText = undefined;
		if (this.index < this.overlays.length) {
			overlayText = this.overlays[this.index];
			if (!overlayText) throw 'here';
		} else {
			if (this.overlays.length > this.maxOverlays) return;
			var div = document.createElement('div');
			div.style.position = 'absolute';
			div.style.pointerEvents = 'none';
			div.style.zIndex = 1;
			overlayText = {div:div};
			this.overlays.push(overlayText); 
		}

		var pos = [];
		vec3.sub(pos, target.pos, orbitTarget.pos);
		pos[3] = 1;
		vec4.transformMat4(pos, pos, glutil.scene.mvMat);
var distInM = vec3.length(pos);
		vec4.transformMat4(pos, pos, glutil.scene.projMat);
		vec4.scale(pos, pos, 1/pos[3]);
		if (pos[0] < -1 || pos[0] > 1 || pos[1] < -1 || pos[1] > 1 || pos[2] < -1 || pos[2] > 1) return;
		
		var sx = parseInt((1+pos[0])/2 * canvas.width);
		var sy = parseInt((1-pos[1])/2 * canvas.height);
		//if (sx == overlayText.x && sy == overlayText.y) return;

var distInParsecs = distInM / metersPerUnits.pc;
var apparentMagnitude = target.magnitude + 5 * (Math.log10(distInParsecs) - 1)
		
		$(overlayText.div).text(
			target.name	//+' '+apparentMagnitude.toFixed(4)
		);

		overlayText.target = target;
		overlayText.x = sx;
		overlayText.y = sy;
		overlayText.div.style.left = sx;
		overlayText.div.style.top = sy;
		document.body.appendChild(overlayText.div);
		var sw = overlayText.div.clientWidth;
		var sh = overlayText.div.clientHeight;
		
		//bump out of the way of other texts
		var bumpedEver = true;
		sx -= sw/2;
		sy -= sh/2;

		var startSX = sx;
		var startSY = sy;
		for (var tries = 0; tries < this.index; ++tries) {
			var bumped = false;
			for (var i = 0; i <  this.index; ++i) {
				var o = this.overlays[i];
				var overlapX = sx < o.x + o.div.clientWidth && sx + sw > o.x;
				var overlapY = sy < o.y + o.div.clientHeight && sy + sh > o.y;
				if (overlapX && overlapY) {
					var push = 0;	//push direction.  0 = x, 1 = y 
					/*push by distance from collided object*/
					var dx = (sx + .5 * sw) - (o.x + .5 * o.div.clientWidth);
					var dy = (sy + .5 * sh) - (o.x + .5 * o.div.clientHeight);
					/**/
					/*push by distance from start * /
					var dx = (sx + .5 * sw) - startSX;
					var dy = (sy + .5 * sh) - startSY;
					/**/
					if (overlapX && overlapY) {
						if (Math.abs(dy) > Math.abs(dx)) push = 1;
					} else if (overlapY) {
						push = 1;
					}
					var padding = 1;
					if (push == 0) {
						if (dx > 0) {
							sx = o.x + o.div.clientWidth + padding;
						} else {
							sx = o.x - sw - padding;
						}
						bumped = true;
						bumpedEver = true;
					} else {	//push == 1
						if (dy > 0) {
							sy = o.y + o.div.clientHeight + padding;
						} else {
							sy = o.y - sh - padding;
						}
						bumped = true;
						bumpedEver = true;
					}
				}
			}
			if (!bumped) break;
		}
		if (bumpedEver) {
			overlayText.x = sx;
			overlayText.y = sy;
			overlayText.div.style.left = sx;
			overlayText.div.style.top = sy;
		}
		
		++this.index;
	};
};
