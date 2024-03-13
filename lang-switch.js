// language-switch.js

// 获取浏览器语言
const userLang = navigator.language || navigator.userLanguage;

// 根据语言选择加载对应的页面
if (userLang.startsWith('zh')) {
  // 如果是中文，跳转到中文版页面
  window.location.href = 'index-zh.html';
} else {
  // 否则跳转到英文版页面
  window.location.href = 'index.html';
}
