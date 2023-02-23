(function(){function r(e,n,t){function o(i,f){if(!n[i]){if(!e[i]){var c="function"==typeof require&&require;if(!f&&c)return c(i,!0);if(u)return u(i,!0);var a=new Error("Cannot find module '"+i+"'");throw a.code="MODULE_NOT_FOUND",a}var p=n[i]={exports:{}};e[i][0].call(p.exports,function(r){var n=e[i][1][r];return o(n||r)},p,p.exports,r,e,n,t)}return n[i].exports}for(var u="function"==typeof require&&require,i=0;i<t.length;i++)o(t[i]);return o}return r})()({1:[function(require,module,exports){
/* global _satellite */

'use strict';

const webvitals = require('web-vitals');

window.webVitalsData = {};

const View = {
	init: function () {
		const view = this;
		webvitals.onCLS(view.addCLSToQueue);
		webvitals.onFID(view.addFIDToQueue);
		webvitals.onLCP(view.addLCPToQueue);
		webvitals.onFCP(view.addFCPToQueue);
		webvitals.onTTFB(view.addTTFBToQueue);
		webvitals.onINP(view.addINPToQueue);

		addEventListener('error', (event) => {
			view.handleJSError(view, event);
		});

		addEventListener('visibilitychange', () => {
			if (document.visibilityState === 'hidden') {
				view.sendPerformanceData(view);
			}
		});

		// support for Safari
		addEventListener('pagehide', () => {
			view.sendPerformanceData(view);
		});
	},
	addCLSToQueue: function (cls) {
		window.webVitalsData.cls = cls.value.toFixed(3);
		// only tracking the first one, once that gets fixed then if there is a second one that will start showing up
		if (cls.value > 0) {
			const largestCLS = cls.entries.reduce((x, y) => x.value > y.value ? x : y);
			if (largestCLS.sources.length > 0) {
				window.webVitalsData.clsSourceNode = largestCLS.sources[0].node;
			}
		}
	},
	addFIDToQueue: function (fid) {
		window.webVitalsData.fid = Math.round(fid.value);
		// only the first input matters for first input delay so no need to look beyond first item in array
		window.webVitalsData.fidAction = fid.entries[0].name;
		window.webVitalsData.fidTargetNode = fid.entries[0].target;
	},
	addLCPToQueue: function (lcp) {
		window.webVitalsData.lcp = Math.round(lcp.value);
		// the last item in the array is always the largest, the other items were the largest until the next one came along
		if (lcp.entries.length) {
			const lcpEntry = lcp.entries[lcp.entries.length - 1];
			window.webVitalsData.lcpUrl = lcpEntry.url;
			window.webVitalsData.lcpTargetNode = lcpEntry.element;
		}
	},
	addFCPToQueue: function (fcp) {
		window.webVitalsData.fcp = Math.round(fcp.value);
	},
	addTTFBToQueue: function (ttfb) {
		window.webVitalsData.ttfb = Math.round(ttfb.value);
		if (ttfb.entries.length > 0) {
			window.webVitalsData.domComplete = Math.round(ttfb.entries[0].domComplete);
			window.webVitalsData.domInteractive = Math.round(ttfb.entries[0].domInteractive);
			window.webVitalsData.entryType = ttfb.entries[0].entryType;
		}
	},
	addINPToQueue: function (inp) {
		window.webVitalsData.inp = inp.value;
		if (inp.entries.length) {
			// picking the largest, all are issues but fixing the largest first will go a long way
			const event = inp.entries.reduce((x, y) => x.duration > y.duration ? x : y, {duration: 0});
			window.webVitalsData.inpName = event.name;
			window.webVitalsData.inpDuration = event.duration;
			if (event.target) {
				window.webVitalsData.inpTargetNode = event.target;
			}
		}
	},
	htmlElementToString: function (element) {
		if (element) {
			let stringName = element.localName;
			if (element.id && element.id.length > 0) {
				stringName += '#' + element.id;
			}
			if (element.className && element.className.length > 0) {
				stringName += '.' + element.className.trim().replaceAll(/\s+/g, '.');
			}
			return stringName;
		}
		return '';
	},
	htmlElementToAncestorString: function (view, element) {
		let ancestorString = '';
		while (element && element?.localName !== 'body') {
			if (ancestorString) {
				ancestorString = ' > ' + ancestorString;
			}
			ancestorString = view.htmlElementToString(element) + ancestorString;
			element = element.parentElement;
		}
		return ancestorString;
	},
	generateSessionKey() {
		if (typeof(sessionStorage.cwvSessionId) === 'undefined') {
			sessionStorage.cwvSessionId = 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, function (c) {
				const r = Math.random() * 16 | 0,
					v = c === 'x' ? r : r & 0x3 | 0x8;
				return v.toString(16);
			});
		}
		return sessionStorage.cwvSessionId;
	},
	flushData() {
		window.webVitalsData = {};
	},
	sendPerformanceData(view) {
		if (Object.keys(window.webVitalsData).length && !location.hostname.includes('aem')) {
			if (Object.hasOwn) {
				if (Object.hasOwn(window.webVitalsData, 'clsSourceNode')) {
					window.webVitalsData.clsSource = view.htmlElementToString(window.webVitalsData.clsSourceNode);
					window.webVitalsData.clsSourceTree = view.htmlElementToAncestorString(view, window.webVitalsData.clsSourceNode);
					delete window.webVitalsData.clsSourceNode;
				}
				if (Object.hasOwn(window.webVitalsData, 'fidTargetNode')) {
					window.webVitalsData.fidTarget = view.htmlElementToString(window.webVitalsData.fidTargetNode);
					window.webVitalsData.fidTargetTree = view.htmlElementToAncestorString(view, window.webVitalsData.fidTargetNode);
					delete window.webVitalsData.fidTargetNode;
				}
				if (Object.hasOwn(window.webVitalsData, 'lcpTargetNode')) {
					window.webVitalsData.lcpTarget = view.htmlElementToString(window.webVitalsData.lcpTargetNode);
					window.webVitalsData.lcpTargetTree = view.htmlElementToAncestorString(view, window.webVitalsData.lcpTargetNode);
					delete window.webVitalsData.lcpTargetNode;
				}
				if (Object.hasOwn(window.webVitalsData, 'inpTargetNode')) {
					window.webVitalsData.inpTarget = view.htmlElementToString(window.webVitalsData.inpTargetNode);
					window.webVitalsData.inpTargetTree = view.htmlElementToAncestorString(view, window.webVitalsData.inpTargetNode);
					delete window.webVitalsData.inpTargetNode;
				}
			}
			view.sendData(window.webVitalsData, 'web vitals');
			view.flushData();
		}
	},
	handleJSError(view, event) {
		let errorData = {};

		errorData.errorMessage = event.message;
		if (event.error) {
			errorData.errorName = event.error.name;
		}

		try { // non-standard features but supported in all browsers
			errorData.errorFilename = event.filename;
			errorData.lineNumber = event.lineno;
			errorData.columnNumber = event.colno;
			if (event.error) {
				errorData.errorStack = event.error.stack;
				errorData.errorMessage = event.error.message;
			}
		} catch (e) { }
		view.sendData(errorData, 'javascript');
	},
	sendData(data, dataType) {
		let dataKey = 'WEB_CORE_WEB_VITALS';
		if (dataType === 'javascript') {
			dataKey = 'JS_ERRORS';
		}
		data.country = window.MathWorksCountryCode;
		data.region = window.MathWorksRegionCode;
		data.userAgent = navigator.userAgent;
		data.query = location.search;
		data.host = location.host;
		data.pathname = location.pathname;
		if (typeof _satellite !== 'undefined') {
			data.sessionid = _satellite.getVisitorId().getMarketingCloudVisitorID();
		}

		let url = 'https://udc-service.mathworks.com/udc/service/v1/events';
		if (location.hostname.includes('-')) {
			url = 'https://udc-service-integ3.mathworks.com/udc/service/v1/events';
		}
		fetch(url, {
			method: 'POST',
			headers: {
				'Content-Type': 'application/json',
				'Accept': 'application/json',
				'x-mw-udc-client-version': '1.0',
				'x-mw-udc-application-name': 'WEB_PAGE_METRICS',
				'x-mw-udc-application-version': '0.5',
				'x-mw-authentication': '2369b91c-a64f-11ec-b909-0242ac120002'
			},
			cache: 'no-cache',
			keepalive: true,
			body: JSON.stringify({
				'Event': [{
					'sessionKey': this.generateSessionKey(),
					'eventKey': dataKey,
					'eventDate': new Date().toISOString().slice(0, 23),
					'eventData': JSON.stringify(data)
				}]
			})
		});
	}
};

module.exports = View;
View.init();
},{"web-vitals":2}],2:[function(require,module,exports){
!function(e,n){"object"==typeof exports&&"undefined"!=typeof module?n(exports):"function"==typeof define&&define.amd?define(["exports"],n):n((e="undefined"!=typeof globalThis?globalThis:e||self).webVitals={})}(this,(function(e){"use strict";var n,t,r,i,o,a=-1,c=function(e){addEventListener("pageshow",(function(n){n.persisted&&(a=n.timeStamp,e(n))}),!0)},u=function(){return window.performance&&performance.getEntriesByType&&performance.getEntriesByType("navigation")[0]},f=function(){var e=u();return e&&e.activationStart||0},s=function(e,n){var t=u(),r="navigate";return a>=0?r="back-forward-cache":t&&(r=document.prerendering||f()>0?"prerender":t.type.replace(/_/g,"-")),{name:e,value:void 0===n?-1:n,rating:"good",delta:0,entries:[],id:"v3-".concat(Date.now(),"-").concat(Math.floor(8999999999999*Math.random())+1e12),navigationType:r}},d=function(e,n,t){try{if(PerformanceObserver.supportedEntryTypes.includes(e)){var r=new PerformanceObserver((function(e){n(e.getEntries())}));return r.observe(Object.assign({type:e,buffered:!0},t||{})),r}}catch(e){}},l=function(e,n){var t=function t(r){"pagehide"!==r.type&&"hidden"!==document.visibilityState||(e(r),n&&(removeEventListener("visibilitychange",t,!0),removeEventListener("pagehide",t,!0)))};addEventListener("visibilitychange",t,!0),addEventListener("pagehide",t,!0)},p=function(e,n,t,r){var i,o;return function(a){n.value>=0&&(a||r)&&((o=n.value-(i||0))||void 0===i)&&(i=n.value,n.delta=o,n.rating=function(e,n){return e>n[1]?"poor":e>n[0]?"needs-improvement":"good"}(n.value,t),e(n))}},v=-1,m=function(){return"hidden"!==document.visibilityState||document.prerendering?1/0:0},h=function(){l((function(e){var n=e.timeStamp;v=n}),!0)},g=function(){return v<0&&(v=m(),h(),c((function(){setTimeout((function(){v=m(),h()}),0)}))),{get firstHiddenTime(){return v}}},y=function(e,n){n=n||{};var t,r=[1800,3e3],i=g(),o=s("FCP"),a=function(e){e.forEach((function(e){"first-contentful-paint"===e.name&&(l&&l.disconnect(),e.startTime<i.firstHiddenTime&&(o.value=e.startTime-f(),o.entries.push(e),t(!0)))}))},u=window.performance&&window.performance.getEntriesByName&&window.performance.getEntriesByName("first-contentful-paint")[0],l=u?null:d("paint",a);(u||l)&&(t=p(e,o,r,n.reportAllChanges),u&&a([u]),c((function(i){o=s("FCP"),t=p(e,o,r,n.reportAllChanges),requestAnimationFrame((function(){requestAnimationFrame((function(){o.value=performance.now()-i.timeStamp,t(!0)}))}))})))},T=!1,E=-1,C=function(e,n){n=n||{};var t=[.1,.25];T||(y((function(e){E=e.value})),T=!0);var r,i=function(n){E>-1&&e(n)},o=s("CLS",0),a=0,u=[],f=function(e){e.forEach((function(e){if(!e.hadRecentInput){var n=u[0],t=u[u.length-1];a&&e.startTime-t.startTime<1e3&&e.startTime-n.startTime<5e3?(a+=e.value,u.push(e)):(a=e.value,u=[e]),a>o.value&&(o.value=a,o.entries=u,r())}}))},v=d("layout-shift",f);v&&(r=p(i,o,t,n.reportAllChanges),l((function(){f(v.takeRecords()),r(!0)})),c((function(){a=0,E=-1,o=s("CLS",0),r=p(i,o,t,n.reportAllChanges)})))},w={passive:!0,capture:!0},L=new Date,b=function(e,i){n||(n=i,t=e,r=new Date,F(removeEventListener),S())},S=function(){if(t>=0&&t<r-L){var e={entryType:"first-input",name:n.type,target:n.target,cancelable:n.cancelable,startTime:n.timeStamp,processingStart:n.timeStamp+t};i.forEach((function(n){n(e)})),i=[]}},A=function(e){if(e.cancelable){var n=(e.timeStamp>1e12?new Date:performance.now())-e.timeStamp;"pointerdown"==e.type?function(e,n){var t=function(){b(e,n),i()},r=function(){i()},i=function(){removeEventListener("pointerup",t,w),removeEventListener("pointercancel",r,w)};addEventListener("pointerup",t,w),addEventListener("pointercancel",r,w)}(n,e):b(n,e)}},F=function(e){["mousedown","keydown","touchstart","pointerdown"].forEach((function(n){return e(n,A,w)}))},I=function(e,r){r=r||{};var o,a=[100,300],u=g(),f=s("FID"),v=function(e){e.startTime<u.firstHiddenTime&&(f.value=e.processingStart-e.startTime,f.entries.push(e),o(!0))},m=function(e){e.forEach(v)},h=d("first-input",m);o=p(e,f,a,r.reportAllChanges),h&&l((function(){m(h.takeRecords()),h.disconnect()}),!0),h&&c((function(){var c;f=s("FID"),o=p(e,f,a,r.reportAllChanges),i=[],t=-1,n=null,F(addEventListener),c=v,i.push(c),S()}))},P=0,M=1/0,k=0,B=function(e){e.forEach((function(e){e.interactionId&&(M=Math.min(M,e.interactionId),k=Math.max(k,e.interactionId),P=k?(k-M)/7+1:0)}))},D=function(){return o?P:performance.interactionCount||0},x=function(){"interactionCount"in performance||o||(o=d("event",B,{type:"event",buffered:!0,durationThreshold:0}))},N=0,R=function(){return D()-N},q=[],H={},O=function(e){var n=q[q.length-1],t=H[e.interactionId];if(t||q.length<10||e.duration>n.latency){if(t)t.entries.push(e),t.latency=Math.max(t.latency,e.duration);else{var r={id:e.interactionId,latency:e.duration,entries:[e]};H[r.id]=r,q.push(r)}q.sort((function(e,n){return n.latency-e.latency})),q.splice(10).forEach((function(e){delete H[e.id]}))}},j=function(e,n){n=n||{};var t=[200,500];x();var r,i=s("INP"),o=function(e){e.forEach((function(e){(e.interactionId&&O(e),"first-input"===e.entryType)&&(!q.some((function(n){return n.entries.some((function(n){return e.duration===n.duration&&e.startTime===n.startTime}))}))&&O(e))}));var n,t=(n=Math.min(q.length-1,Math.floor(R()/50)),q[n]);t&&t.latency!==i.value&&(i.value=t.latency,i.entries=t.entries,r())},a=d("event",o,{durationThreshold:n.durationThreshold||40});r=p(e,i,t,n.reportAllChanges),a&&(a.observe({type:"first-input",buffered:!0}),l((function(){o(a.takeRecords()),i.value<0&&R()>0&&(i.value=0,i.entries=[]),r(!0)})),c((function(){q=[],N=D(),i=s("INP"),r=p(e,i,t,n.reportAllChanges)})))},_={},V=function(e,n){n=n||{};var t,r=[2500,4e3],i=g(),o=s("LCP"),a=function(e){var n=e[e.length-1];if(n){var r=n.startTime-f();r<i.firstHiddenTime&&(o.value=r,o.entries=[n],t())}},u=d("largest-contentful-paint",a);if(u){t=p(e,o,r,n.reportAllChanges);var v=function(){_[o.id]||(a(u.takeRecords()),u.disconnect(),_[o.id]=!0,t(!0))};["keydown","click"].forEach((function(e){addEventListener(e,v,{once:!0,capture:!0})})),l(v,!0),c((function(i){o=s("LCP"),t=p(e,o,r,n.reportAllChanges),requestAnimationFrame((function(){requestAnimationFrame((function(){o.value=performance.now()-i.timeStamp,_[o.id]=!0,t(!0)}))}))}))}},z=function e(n){document.prerendering?addEventListener("prerenderingchange",(function(){return e(n)}),!0):"complete"!==document.readyState?addEventListener("load",(function(){return e(n)}),!0):setTimeout(n,0)},G=function(e,n){n=n||{};var t=[800,1800],r=s("TTFB"),i=p(e,r,t,n.reportAllChanges);z((function(){var o=u();if(o){if(r.value=Math.max(o.responseStart-f(),0),r.value<0||r.value>performance.now())return;r.entries=[o],i(!0),c((function(){r=s("TTFB",0),(i=p(e,r,t,n.reportAllChanges))(!0)}))}}))};e.getCLS=C,e.getFCP=y,e.getFID=I,e.getINP=j,e.getLCP=V,e.getTTFB=G,e.onCLS=C,e.onFCP=y,e.onFID=I,e.onINP=j,e.onLCP=V,e.onTTFB=G,Object.defineProperty(e,"__esModule",{value:!0})}));

},{}]},{},[1]);
