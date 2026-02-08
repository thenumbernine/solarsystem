import {vec3, vec4} from '/js/gl-matrix-3.4.1/index.js';
import {cfg} from './globals.js';
import {ui} from './ui.js';
import {metersPerUnits} from './units.js';

class OverlayTexts {
	constructor() {
		this.overlays = [];
		this.index = 0;
		this.maxOverlays = 20;
	}

	updateBegin() {
		this.index = 0;
	}

	updateEnd() {
		for (let i = this.index; i < this.overlays.length; ++i) {
			const div = this.overlays[i].div;
			const parent = div.parentNode;
			if (parent) parent.removeChild(div);
		}
	}

	add(target) {
		if (!cfg.showNames) return;

		for (let i = 0; i < this.index; ++i) {
			if (this.overlays[i].target == target) return;
		}

		let overlayText = undefined;
		if (this.index < this.overlays.length) {
			overlayText = this.overlays[this.index];
			if (!overlayText) throw 'here';
		} else {
			if (this.overlays.length > this.maxOverlays) return;
			let div = document.createElement('div');
			div.style.position = 'absolute';
			div.style.pointerEvents = 'none';
			div.style.zIndex = 1;
			overlayText = {div:div};
			this.overlays.push(overlayText);
		}

		let pos = [];
		vec3.sub(pos, target.pos, cfg.orbitTarget.pos);
		pos[3] = 1;
		vec4.transformMat4(pos, pos, glutil.scene.mvMat);
let distInM = vec3.length(pos);
		vec4.transformMat4(pos, pos, glutil.scene.projMat);
		vec4.scale(pos, pos, 1/pos[3]);
		if (pos[0] < -1 || pos[0] > 1 || pos[1] < -1 || pos[1] > 1 || pos[2] < -1 || pos[2] > 1) return;

		const canvas = ui.canvas;
		let sx = parseInt((1+pos[0])/2 * canvas.width);
		let sy = parseInt((1-pos[1])/2 * canvas.height);
		//if (sx == overlayText.x && sy == overlayText.y) return;

let distInParsecs = distInM / metersPerUnits.pc;
let apparentMagnitude = target.magnitude + 5 * (Math.log10(distInParsecs) - 1)

		overlayText.div.innerText =
			target.name	//+' '+apparentMagnitude.toFixed(4)
		;

		overlayText.target = target;
		overlayText.x = sx;
		overlayText.y = sy;
		overlayText.div.style.left = sx+'px';
		overlayText.div.style.top = sy+'px';
		document.body.appendChild(overlayText.div);
		let sw = overlayText.div.clientWidth;
		let sh = overlayText.div.clientHeight;

		//bump out of the way of other texts
		let bumpedEver = true;
		sx -= sw/2;
		sy -= sh/2;

		let startSX = sx;
		let startSY = sy;
		for (let tries = 0; tries < this.index; ++tries) {
			let bumped = false;
			for (let i = 0; i <  this.index; ++i) {
				let o = this.overlays[i];
				let overlapX = sx < o.x + o.div.clientWidth && sx + sw > o.x;
				let overlapY = sy < o.y + o.div.clientHeight && sy + sh > o.y;
				if (overlapX && overlapY) {
					let push = 0;	//push direction.  0 = x, 1 = y
					/*push by distance from collided object*/
					let dx = (sx + .5 * sw) - (o.x + .5 * o.div.clientWidth);
					let dy = (sy + .5 * sh) - (o.x + .5 * o.div.clientHeight);
					/**/
					/*push by distance from start * /
					let dx = (sx + .5 * sw) - startSX;
					let dy = (sy + .5 * sh) - startSY;
					/**/
					if (overlapX && overlapY) {
						if (Math.abs(dy) > Math.abs(dx)) push = 1;
					} else if (overlapY) {
						push = 1;
					}
					let padding = 1;
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
			overlayText.div.style.left = sx+'px';
			overlayText.div.style.top = sy+'px';
		}

		++this.index;
	}
}
const overlayTexts = new OverlayTexts();
export {overlayTexts};
