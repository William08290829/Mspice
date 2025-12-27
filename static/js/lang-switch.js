// language-switch.js

document.addEventListener('DOMContentLoaded', () => {
  const userLang = navigator.language || navigator.userLanguage;
  const prefLang = localStorage.getItem('user_pref_lang');

  // Only redirect if no preference is saved
  if (!prefLang) {
    if (userLang.startsWith('zh')) {
      window.location.href = 'index-zh.html';
    } else {
      // Already on index.html, no action needed for EN
      // If we were on index-zh.html, we'd redirect to index.html here
    }
  }

  // Attach click listeners to language options to save preference
  const langOptions = document.querySelectorAll('.lang-option');
  langOptions.forEach(option => {
    option.addEventListener('click', (e) => {
      // Determine language based on link destination or text
      // Here 'index-zh.html' implies Chinese
      if (option.getAttribute('href').includes('zh')) {
        localStorage.setItem('user_pref_lang', 'zh');
      } else {
        localStorage.setItem('user_pref_lang', 'en');
      }
    });
  });
});
